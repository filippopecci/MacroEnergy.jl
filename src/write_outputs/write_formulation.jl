###############################################################################
# write_formulation.jl
#
# Generates a compact, human-readable mathematical formulation of the model that
# is actually built for a given `System`/`Case`, in both LaTeX (`.tex`, -> PDF)
# and Markdown (`.md`, Documenter/KaTeX) form.
#
# Math content comes from a dedicated, dispatch-based API:
#   - `formulation(constraint, element)` for constraints stored as
#     `AbstractTypeConstraint` instances (mirrors `add_model_constraint!`);
#   - `structural_formulations(element)` for relations added directly with
#     `@constraint(...)` during model build;
#   - balance constraints are *reconstructed* from `balance_data` so the actual
#     signed flow terms (and demand offsets) are shown, not an opaque expression.
#
# A compact symbol notation (see `formulation_parameters`) keeps equations short
# enough to fit the page; symbols are listed in a Notation section of the output.
###############################################################################

# -----------------------------------------------------------------------------
# Data structures
# -----------------------------------------------------------------------------

"""
    ConstraintFormulation

A render-target-agnostic description of a constraint's math. `body` holds one or
more LaTeX equation strings (no fences); each is rendered in an `aligned`/`align`
environment by the writers.
"""
Base.@kwdef struct ConstraintFormulation
    name::String
    body::Vector{String} = String[]
    indexsets::Vector{String} = String[]
    symbols::Vector{Pair{String,String}} = Pair{String,String}[]
    conditions::Vector{String} = String[]
end

"""
    VariableFormulation

A structured description of a single decision-variable family.
"""
Base.@kwdef struct VariableFormulation
    symbol::String
    name::String
    domain::String = ""
    indexset::String = ""
    condition::String = ""
end

"""
    ActiveConstraintEntry

A `ConstraintFormulation` together with the element ids it applies to and the
planning periods in which it is active. Used to print each constraint once
(compactly) with an "applies to" list.
"""
Base.@kwdef mutable struct ActiveConstraintEntry
    formulation::ConstraintFormulation
    element_kind::String
    category::Int                 # 1 Planning, 2 Operation, 3 Policy, 4 Structural
    elements::Vector{Symbol} = Symbol[]
    periods::Vector{Int} = Int[]
end

# -----------------------------------------------------------------------------
# Compact notation: parameters (variables are listed by `formulation_variables`)
# -----------------------------------------------------------------------------

"""
    formulation_parameters()

Symbol => meaning for the parameters used in the formulation bodies. Listed in
the Notation section of the output so the short symbols are legible.
"""
formulation_parameters() = Pair{String,String}[
    raw"\alpha_{e,t}"                => "availability (capacity factor) of edge e at time t",
    raw"K_{e}"                       => "capacity per built unit of edge e",
    raw"x^{\mathrm{ex}}_{y}"         => "existing capacity of edge/storage y",
    raw"\underline{x}_{y},\ \overline{x}_{y}" => "minimum / maximum capacity of y",
    raw"\overline{x}^{\,\mathrm{new}}_{y}"    => "maximum new capacity of y",
    raw"d_{n,t}"                     => "demand at node n, time t",
    raw"\varphi_{e}"                 => "minimum flow fraction of edge e",
    raw"r^{\uparrow}_{e},\ r^{\downarrow}_{e}" => "ramp up / down fraction of edge e",
    raw"T^{\uparrow}_{e},\ T^{\downarrow}_{e}" => "minimum up / down time of edge e",
    raw"\eta_{e}"                    => "efficiency of edge e",
    raw"\underline{\sigma}_{g},\ \overline{\sigma}_{g}" => "min / max storage level (fraction of capacity) of g",
    raw"\underline{\delta}_{g},\ \overline{\delta}_{g}" => "min / max storage duration of g",
    raw"\gamma_{g}"                  => "charge/discharge capacity ratio of storage g",
    raw"\lambda_{g}"                 => "self-discharge (loss) fraction of storage g",
    raw"\varphi^{\mathrm{out}}_{g}"  => "minimum outflow fraction of storage g",
    raw"\mathrm{reg}_{e,t},\ \mathrm{rsv}_{e,t}" => "regulation / reserve terms on edge e",
    raw"b_{n}"                       => "policy right-hand side at node n",
    raw"\varepsilon_{n,t}"           => "net emissions balance at node n, time t",
    raw"\omega_{t}"                  => "weight of the subperiod containing time t",
    raw"\theta_{r}"                  => "retrofit efficiency of retrofit option r",
]

# -----------------------------------------------------------------------------
# Per-constraint formulation API (dispatched identically to add_model_constraint!)
# -----------------------------------------------------------------------------

# Fallback keeps the writer robust if a constraint has no authored formulation.
function formulation(ct::AbstractTypeConstraint, element)
    return ConstraintFormulation(
        name = string(nameof(typeof(ct))),
        conditions = ["Formulation not yet authored for $(nameof(typeof(ct))) on $(nameof(typeof(element)))."],
    )
end

## --- Planning constraints -----------------------------------------------------

formulation(::MinCapacityConstraint, ::Union{AbstractEdge,AbstractStorage}) = ConstraintFormulation(
    name = "Minimum capacity", body = [raw"x_{y} \geq \underline{x}_{y}"])

formulation(::MaxCapacityConstraint, ::Union{AbstractEdge,AbstractStorage}) = ConstraintFormulation(
    name = "Maximum capacity", body = [raw"x_{y} \leq \overline{x}_{y}"])

formulation(::MaxNewCapacityConstraint, ::Union{AbstractEdge,AbstractStorage}) = ConstraintFormulation(
    name = "Maximum new capacity", body = [raw"x^{\mathrm{new}}_{y} \leq \overline{x}^{\,\mathrm{new}}_{y}"])

formulation(::AgeBasedRetirementConstraint, ::Union{AbstractEdge,AbstractStorage}) = ConstraintFormulation(
    name = "Age-based retirement",
    body = [raw"\sum_{k \leq K^{\mathrm{ret}}(y)} x^{\mathrm{new}}_{y,k} + \underline{x}^{\mathrm{ret}}_{y} \;\leq\; \sum_{k} \big( x^{\mathrm{ret}}_{y,k} + x^{\mathrm{rf}}_{y,k} \big)"],
    symbols = [raw"K^{\mathrm{ret}}(y)" => "retirement period of y", raw"x^{\mathrm{new}}_{y,k}" => "capacity of y built in period k"],
    conditions = ["all capacity built up to the retirement period must retire or be retrofitted"])

## --- Operation constraints ----------------------------------------------------

formulation(::CapacityConstraint, ::UnidirectionalEdge) = ConstraintFormulation(
    name = "Capacity (unidirectional edge)",
    body = [raw"f_{e,t} \leq \alpha_{e,t}\, x_{e}"],
    indexsets = [raw"e \in \mathcal{E},\ t \in \mathcal{T}"], conditions = ["only if e has capacity"])

formulation(::CapacityConstraint, ::BidirectionalEdge) = ConstraintFormulation(
    name = "Capacity (bidirectional edge)",
    body = [raw"i\, f_{e,t} \leq \alpha_{e,t}\, x_{e}"],
    indexsets = [raw"i \in \{-1, 1\},\ t \in \mathcal{T}"], conditions = ["only if e has capacity"])

formulation(::CapacityConstraint, ::EdgeWithUC) = ConstraintFormulation(
    name = "Capacity (edge with unit commitment)",
    body = [raw"f_{e,t} \leq \alpha_{e,t}\, K_{e}\, u_{e,t}"],
    indexsets = [raw"t \in \mathcal{T}"])

formulation(::MustRunConstraint, ::AbstractEdge) = ConstraintFormulation(
    name = "Must run", body = [raw"f_{e,t} = \alpha_{e,t}\, x_{e}"],
    indexsets = [raw"t \in \mathcal{T}"], conditions = ["unidirectional edges with capacity only"])

formulation(::MinFlowConstraint, ::UnidirectionalEdge) = ConstraintFormulation(
    name = "Minimum flow", body = [raw"f_{e,t} \geq \varphi_{e}\, x_{e}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::MinFlowConstraint, ::EdgeWithUC) = ConstraintFormulation(
    name = "Minimum flow (edge with unit commitment)",
    body = [raw"f_{e,t} \geq \varphi_{e}\, K_{e}\, u_{e,t}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::RampingLimitConstraint, ::EdgeWithoutUC) = ConstraintFormulation(
    name = "Ramping limits",
    body = [
        raw"f_{e,t} - f_{e,t-1} + \mathrm{reg}_{e,t} + \mathrm{rsv}_{e,t} \;\leq\; r^{\uparrow}_{e}\, x_{e}",
        raw"f_{e,t-1} - f_{e,t} + \mathrm{reg}_{e,t} + \mathrm{rsv}_{e,t} \;\leq\; r^{\downarrow}_{e}\, x_{e}",
    ],
    indexsets = [raw"t \in \mathcal{T}"])

formulation(::RampingLimitConstraint, ::EdgeWithUC) = ConstraintFormulation(
    name = "Ramping limits (edge with unit commitment)",
    body = [
        raw"f_{e,t} - f_{e,t-1} + \mathrm{reg}_{e,t} + \mathrm{rsv}_{e,t} \;\leq\; r^{\uparrow}_{e} K_{e} (u_{e,t} - u^{+}_{e,t}) - m^{\uparrow}_{e,t} K_{e} u^{+}_{e,t} + \varphi_{e} K_{e} u^{-}_{e,t}",
        raw"f_{e,t-1} - f_{e,t} + \mathrm{reg}_{e,t} + \mathrm{rsv}_{e,t} \;\leq\; r^{\downarrow}_{e} K_{e} (u_{e,t} - u^{+}_{e,t}) + \varphi_{e} K_{e} u^{+}_{e,t} - m^{\downarrow}_{e,t} K_{e} u^{-}_{e,t}",
    ],
    indexsets = [raw"t \in \mathcal{T}"],
    symbols = [
        raw"m^{\uparrow}_{e,t}"   => raw"\min(\alpha_{e,t},\ \max(\varphi_{e},\ r^{\uparrow}_{e}))",
        raw"m^{\downarrow}_{e,t}" => raw"\min(\alpha_{e,t},\ \max(\varphi_{e},\ r^{\downarrow}_{e}))",
    ])

formulation(::MinUpTimeConstraint, ::EdgeWithUC) = ConstraintFormulation(
    name = "Minimum up time", body = [raw"u_{e,t} \geq \sum_{h=0}^{T^{\uparrow}_{e}-1} u^{+}_{e,t-h}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::MinDownTimeConstraint, ::EdgeWithUC) = ConstraintFormulation(
    name = "Minimum down time", body = [raw"\frac{x_{e}}{K_{e}} - u_{e,t} \geq \sum_{h=0}^{T^{\downarrow}_{e}-1} u^{-}_{e,t-h}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::StorageCapacityConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Storage capacity", body = [raw"s_{g,t} \leq x_{g}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::StorageSymmetricCapacityConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Storage symmetric capacity",
    body = [raw"f_{d,t} + f_{c,t} \leq x_{d}"], indexsets = [raw"t \in \mathcal{T}"],
    symbols = [raw"d,\ c" => "discharge and charge edges of storage g", raw"x_{d}" => "capacity of the discharge edge"])

formulation(::StorageChargeDischargeRatioConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Storage charge/discharge ratio", body = [raw"\gamma_{g}\, x_{d} = x_{c}"],
    symbols = [raw"d,\ c" => "discharge and charge edges of storage g"])

formulation(::StorageMaxDurationConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Storage maximum duration", body = [raw"x_{g} \leq \overline{\delta}_{g}\, x_{d}"], conditions = ["only if max duration > 0"])

formulation(::StorageMinDurationConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Storage minimum duration", body = [raw"x_{g} \geq \underline{\delta}_{g}\, x_{d}"], conditions = ["only if min duration > 0"])

formulation(::MaxStorageLevelConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Maximum storage level", body = [raw"s_{g,t} \leq \overline{\sigma}_{g}\, x_{g}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::MinStorageLevelConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Minimum storage level", body = [raw"s_{g,t} \geq \underline{\sigma}_{g}\, x_{g}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::MinStorageOutflowConstraint, ::AbstractStorage) = ConstraintFormulation(
    name = "Minimum storage outflow",
    body = [raw"f_{\mathrm{sp},t} + f_{d,t} \geq \varphi^{\mathrm{out}}_{g}\, x_{d}"], indexsets = [raw"t \in \mathcal{T}"],
    symbols = [raw"\mathrm{sp},\ d" => "spillage and discharge edges of storage g"], conditions = ["HydroRes assets with a spillage edge"])

formulation(::StorageChargeLimitConstraint, ::AbstractEdge) = ConstraintFormulation(
    name = "Storage charge limit",
    body = [raw"\eta_{e}\, f_{e,t} \leq x_{g} - s_{g,t-1}"], indexsets = [raw"t \in \mathcal{T}"],
    symbols = [raw"g" => "storage at the end vertex of e"], conditions = ["edges whose end vertex is a storage"])

formulation(::StorageDischargeLimitConstraint, ::AbstractEdge) = ConstraintFormulation(
    name = "Storage discharge limit",
    body = [raw"\frac{f_{e,t}}{\eta_{e}} \leq s_{g,t-1}"], indexsets = [raw"t \in \mathcal{T}"],
    symbols = [raw"g" => "storage at the start vertex of e"], conditions = ["edges whose start vertex is a storage"])

formulation(::MaxNonServedDemandConstraint, ::Node) = ConstraintFormulation(
    name = "Maximum non-served demand", body = [raw"\sum_{s} \nu_{n,s,t} \leq d_{n,t}"], indexsets = [raw"t \in \mathcal{T}"])

formulation(::MaxNonServedDemandPerSegmentConstraint, ::Node) = ConstraintFormulation(
    name = "Maximum non-served demand per segment",
    body = [raw"\nu_{n,s,t} \leq \overline{\nu}_{n,s}\, d_{n,t}"], indexsets = [raw"s \in \mathcal{S}_n,\ t \in \mathcal{T}"],
    symbols = [raw"\overline{\nu}_{n,s}" => "max non-served demand fraction of segment s", raw"\mathcal{S}_n" => "non-served-demand segments of n"])

formulation(::LongDurationStorageChangeConstraint, ::LongDurationStorage) = ConstraintFormulation(
    name = "Long-duration storage change",
    body = [raw"s^{0}_{g,w} = s_{g,\mathrm{end}(w)} - \Delta_{g,w}"], indexsets = [raw"w \in \mathcal{W}"],
    symbols = [raw"\mathrm{end}(w)" => "last timestep of subperiod w"])

formulation(::LongDurationStorageImplicitMinMaxConstraint, ::LongDurationStorage) = ConstraintFormulation(
    name = "Long-duration storage implicit min/max",
    body = [
        raw"B_{g,p} + S^{\max}_{g,r} - s_{g,\mathrm{t_0}(p)} \leq x_{g}",
        raw"B_{g,p} + S^{\min}_{g,r} - s_{g,\mathrm{t_0}(p)} \geq 0",
    ],
    indexsets = [raw"p \in \text{non-representative subperiods},\ r = \text{rep. subperiod of } p"],
    symbols = [
        raw"B_{g,p}" => raw"(1-\lambda_{g})\, s^{0}_{g,p} + \tfrac{1}{\eta_{d}} f_{d,\mathrm{t_0}(p)} - \eta_{c} f_{c,\mathrm{t_0}(p)}\ \text{(storage balance)}",
        raw"S^{\min}_{g,r},\ S^{\max}_{g,r}" => "min / max storage level over representative subperiod r",
        raw"\mathrm{t_0}(p)" => "first timestep of the representative subperiod modeling p",
    ])

formulation(::MaxInitStorageLevelConstraint, ::LongDurationStorage) = ConstraintFormulation(
    name = "Maximum initial storage level (long-duration)", body = [raw"s^{0}_{g} \leq \overline{\sigma}_{g}\, x_{g}"], conditions = ["Benders decomposition only"])

formulation(::MinInitStorageLevelConstraint, ::LongDurationStorage) = ConstraintFormulation(
    name = "Minimum initial storage level (long-duration)", body = [raw"s^{0}_{g} \geq \underline{\sigma}_{g}\, x_{g}"], conditions = ["Benders decomposition only"])

## --- Policy constraints -------------------------------------------------------

formulation(::CO2CapConstraint, ::Node) = ConstraintFormulation(
    name = "CO₂ cap", body = [raw"\sum_{t \in \mathcal{T}} \varepsilon_{n,t} - \xi_{n} \leq b_{n}"],
    conditions = ["slack ξ present only if an unmet-policy price is set (penalised in the objective)"])

formulation(::CO2StorageConstraint, ::Node) = ConstraintFormulation(
    name = "CO₂ storage budget", body = [raw"\sum_{t \in w} \omega_{t}\, \varepsilon^{\mathrm{st}}_{n,t} \leq B_{n,w}"],
    indexsets = [raw"w \in \mathcal{W}"], symbols = [raw"B_{n,w}" => "CO₂ storage budget for subperiod w"])

formulation(::AggregatedDemandConstraint, ::Node) = ConstraintFormulation(
    name = "Aggregated demand", body = [raw"\sum_{t \in \mathcal{T}} \textstyle\sum_{e \in \mathcal{I}(n)} f_{e,t} \geq b_{n}"],
    symbols = [raw"\mathcal{I}(n)" => "edges flowing into demand node n"])

_constraint_category(c::PolicyConstraint) = 3
_constraint_category(c::OperationConstraint) = 2
_constraint_category(c::PlanningConstraint) = 1
_constraint_category(c::AbstractTypeConstraint) = 2

# -----------------------------------------------------------------------------
# Variable API: `formulation_variables(element)` (compact symbols)
# -----------------------------------------------------------------------------

function formulation_variables(e::AbstractEdge)
    vars = VariableFormulation[]
    push!(vars, VariableFormulation(symbol = raw"f_{e,t}", name = "commodity flow on edge e at time t",
        domain = (e isa BidirectionalEdge ? "free" : raw"\geq 0"), indexset = raw"e \in \mathcal{E},\ t \in \mathcal{T}"))
    if has_capacity(e)
        push!(vars, VariableFormulation(symbol = raw"x_{e}", name = "installed capacity of edge e", domain = raw"\geq 0"))
        can_expand(e) && push!(vars, VariableFormulation(symbol = raw"x^{\mathrm{new}}_{e}", name = "new capacity of edge e", domain = raw"\geq 0"))
        can_retire(e) && push!(vars, VariableFormulation(symbol = raw"x^{\mathrm{ret}}_{e}", name = "retired capacity of edge e", domain = raw"\geq 0"))
        can_retrofit(e) && push!(vars, VariableFormulation(symbol = raw"x^{\mathrm{rf}}_{e}", name = "retrofitted capacity of edge e", domain = raw"\geq 0"))
    end
    if e isa EdgeWithUC
        push!(vars, VariableFormulation(symbol = raw"u_{e,t}", name = "committed units", domain = raw"\geq 0,\ \text{integer}", indexset = raw"e \in \mathcal{E}^{UC},\ t \in \mathcal{T}"))
        push!(vars, VariableFormulation(symbol = raw"u^{+}_{e,t}", name = "started-up units", domain = raw"\geq 0,\ \text{integer}"))
        push!(vars, VariableFormulation(symbol = raw"u^{-}_{e,t}", name = "shut-down units", domain = raw"\geq 0,\ \text{integer}"))
    end
    return vars
end

function formulation_variables(g::AbstractStorage)
    vars = VariableFormulation[]
    push!(vars, VariableFormulation(symbol = raw"s_{g,t}", name = "stored quantity in storage g at time t", domain = raw"\geq 0", indexset = raw"g \in \mathcal{G},\ t \in \mathcal{T}"))
    push!(vars, VariableFormulation(symbol = raw"x_{g}", name = "installed storage capacity", domain = raw"\geq 0"))
    can_expand(g) && push!(vars, VariableFormulation(symbol = raw"x^{\mathrm{new}}_{g}", name = "new storage capacity", domain = raw"\geq 0"))
    can_retire(g) && push!(vars, VariableFormulation(symbol = raw"x^{\mathrm{ret}}_{g}", name = "retired storage capacity", domain = raw"\geq 0"))
    if g isa LongDurationStorage
        push!(vars, VariableFormulation(symbol = raw"s^{0}_{g,w}", name = "initial storage level per subperiod", domain = raw"\geq 0", indexset = raw"g \in \mathcal{G}^{LDS},\ w \in \mathcal{W}"))
        push!(vars, VariableFormulation(symbol = raw"\Delta_{g,w}", name = "storage change over subperiod", domain = "free"))
    end
    return vars
end

function formulation_variables(n::Node)
    vars = VariableFormulation[]
    if !all(max_non_served_demand(n) .== 0)
        push!(vars, VariableFormulation(symbol = raw"\nu_{n,s,t}", name = "non-served demand by segment", domain = raw"\geq 0", indexset = raw"n \in \mathcal{N},\ s,\ t \in \mathcal{T}"))
    end
    if !isempty(supply_segments(n))
        push!(vars, VariableFormulation(symbol = raw"q_{n,s,t}", name = "supply by segment", domain = raw"\geq 0", indexset = raw"n \in \mathcal{N},\ s,\ t \in \mathcal{T}"))
    end
    return vars
end

formulation_variables(::MacroObject) = VariableFormulation[]

# -----------------------------------------------------------------------------
# Structural-constraint API: relations added directly with @constraint(...)
# -----------------------------------------------------------------------------

const _STRUCTURAL_CATEGORY = 4

_unit_count_symbols() = Pair{String,String}[
    raw"x^{\mathrm{new}}_{e}" => raw"K_{e}\, z^{\mathrm{new}}_{e}\ \text{(units built)}",
    raw"x^{\mathrm{ret}}_{e}" => raw"K_{e}\, z^{\mathrm{ret}}_{e}\ \text{(units retired)}",
]

function structural_formulations(e::AbstractEdge)
    fs = ConstraintFormulation[]
    if has_capacity(e)
        if can_retrofit(e)
            push!(fs, ConstraintFormulation(name = "Capacity accounting (retrofittable)",
                body = [raw"x_{e} = x^{\mathrm{new}}_{e} - x^{\mathrm{ret}}_{e} - x^{\mathrm{rf}}_{e} + x^{\mathrm{ex}}_{e}"],
                symbols = _unit_count_symbols(), conditions = ["added automatically for edges with capacity"]))
            push!(fs, ConstraintFormulation(name = "Retirement & retrofit bound",
                body = [raw"x^{\mathrm{rf}}_{e} + x^{\mathrm{ret}}_{e} \leq x^{\mathrm{ex}}_{e}"]))
        else
            push!(fs, ConstraintFormulation(name = "Capacity accounting",
                body = [raw"x_{e} = x^{\mathrm{new}}_{e} - x^{\mathrm{ret}}_{e} + x^{\mathrm{ex}}_{e}"],
                symbols = _unit_count_symbols(), conditions = ["added automatically for edges with capacity"]))
            push!(fs, ConstraintFormulation(name = "Retirement bound", body = [raw"x^{\mathrm{ret}}_{e} \leq x^{\mathrm{ex}}_{e}"]))
        end
    end
    if e isa EdgeWithUC
        push!(fs, ConstraintFormulation(name = "Unit commitment limits",
            body = [raw"u_{e,t},\ u^{+}_{e,t},\ u^{-}_{e,t} \leq x_{e}/K_{e}"], indexsets = [raw"t \in \mathcal{T}"]))
        push!(fs, ConstraintFormulation(name = "Unit commitment state transition",
            body = [raw"u_{e,t} - u_{e,t-1} = u^{+}_{e,t} - u^{-}_{e,t}"], indexsets = [raw"t \in \mathcal{T}"]))
    end
    if e isa BidirectionalEdge && lossy_edge(e)
        push!(fs, ConstraintFormulation(name = "Bidirectional lossy flow split",
            body = [raw"f^{+}_{e,t} - f^{-}_{e,t} = f_{e,t}"], indexsets = [raw"t \in \mathcal{T}"]))
        if has_capacity(e) && any(isa.(e.constraints, CapacityConstraint))
            push!(fs, ConstraintFormulation(name = "Bidirectional lossy capacity",
                body = [raw"f^{+}_{e,t} + f^{-}_{e,t} \leq \alpha_{e,t}\, x_{e}"], indexsets = [raw"t \in \mathcal{T}"]))
        end
    end
    return fs
end

function structural_formulations(g::AbstractStorage)
    fs = ConstraintFormulation[]
    push!(fs, ConstraintFormulation(name = "Storage capacity accounting",
        body = [raw"x_{g} = x^{\mathrm{new}}_{g} - x^{\mathrm{ret}}_{g} + x^{\mathrm{ex}}_{g}"], conditions = ["added automatically"]))
    push!(fs, ConstraintFormulation(name = "Storage retirement bound", body = [raw"x^{\mathrm{ret}}_{g} \leq x^{\mathrm{ex}}_{g}"]))
    if g isa LongDurationStorage
        push!(fs, ConstraintFormulation(name = "Long-duration storage initial-level bound",
            body = [raw"s^{0}_{g,r} \leq x_{g}"], indexsets = [raw"r \in \text{modeled subperiods}"]))
        push!(fs, ConstraintFormulation(name = "Long-duration storage initial-level linking",
            body = [raw"s^{0}_{g,r+1} = s^{0}_{g,r} + \Delta_{g,\,\mathrm{sp}(r)}"], indexsets = [raw"r \in \text{modeled subperiods (cyclic)}"]))
    end
    return fs
end

function structural_formulations(n::Node)
    fs = ConstraintFormulation[]
    if any(isa.(n.constraints, PolicyConstraint))
        push!(fs, ConstraintFormulation(name = "Policy budget balance",
            body = [raw"\sum_{w \in \mathcal{W}} \beta_{n,w} = b_{n}"], symbols = [raw"\beta_{n,w}" => "policy budgeting variable"],
            conditions = ["one per active policy constraint on the node"]))
    end
    if !isempty(supply_segments(n))
        push!(fs, ConstraintFormulation(name = "Supply segment bounds",
            body = [raw"\underline{q}_{n,s,t} \leq q_{n,s,t} \leq \overline{q}_{n,s,t}"], indexsets = [raw"s \in \mathcal{S}^{q}_n,\ t \in \mathcal{T}"],
            conditions = ["each bound applied only where finite (max) or positive (min)"]))
    end
    return fs
end

structural_formulations(::MacroObject) = ConstraintFormulation[]

function system_structural_formulations(system::System)
    out = Tuple{ConstraintFormulation,Vector{Symbol}}[]
    if get(system.settings, :Retrofitting, false)
        retro_edges = [id(e) for e in get_edges(system) if can_retrofit(e) && !ismissing(retrofit_id(e))]
        if !isempty(retro_edges)
            f = ConstraintFormulation(name = "Retrofit capacity linking",
                body = [raw"x^{\mathrm{rf}}_{e} = \sum_{r \in R(e)} \frac{x^{\mathrm{new}}_{r}}{\theta_{r}}"],
                symbols = [raw"R(e)" => "set of retrofit-option edges for edge e"],
                conditions = ["added only when the Retrofitting setting is enabled"])
            push!(out, (f, retro_edges))
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Balance reconstruction: show the actual signed flow terms of each balance
# -----------------------------------------------------------------------------

_balance_op(sense::Symbol) = sense == :eq ? "=" : sense == :le ? raw"\leq" : raw"\geq"

function _edge_incidence(v::AbstractVertex, e::AbstractEdge)
    start_vertex(e) === v && return -1.0   # outgoing
    end_vertex(e) === v && return 1.0      # incoming
    return 1.0
end

# Strip the owning-asset id prefix and "_edge" suffix to get a short edge role.
function _edge_role(asset, e::AbstractEdge)
    s = string(id(e))
    if asset !== nothing
        p = string(id(asset)) * "_"
        startswith(s, p) && (s = s[length(p)+1:end])
    end
    endswith(s, "_edge") && (s = s[1:end-5])
    return s
end

_role_text(role::AbstractString) = "\\text{" * replace(role, "_" => raw"\_") * "}"
_role_flow_sym(role::AbstractString) = "f_{" * _role_text(role) * "}"
# Symbolic conversion coefficient for an edge role (technology-specific value).
_role_coeff_sym(role::AbstractString) = "a_{" * _role_text(role) * "}"

# Reconstruct the signed flow terms of a balance from its BalanceData terms,
# using *symbolic* coefficients (a_role) for any magnitude != 1 so that
# structurally identical balances merge regardless of numeric values. The final
# sign per term is incidence * sign(stored coeff).
# Returns (expression_string, signature_vector, uses_symbolic_coeff::Bool).
function _reconstruct_flow_terms(v::AbstractVertex, asset, bd)
    parts = String[]
    sig = Tuple{String,Int,Bool}[]   # (role, sign, is_unit_magnitude)
    used_symbol = false
    for term in bd.terms
        term.var == :flow || continue
        e = term.obj
        e isa AbstractEdge || continue
        role = _edge_role(asset, e)
        fsym = _role_flow_sym(role)
        inc = _edge_incidence(v, e)
        if term.coeff isa Float64
            coeff = inc * term.coeff
            isunit = isapprox(abs(coeff), 1.0; atol = 1e-9)
            negative = coeff < 0
        else
            isunit = false
            negative = inc < 0
        end
        isunit || (used_symbol = true)
        csym = isunit ? "" : _role_coeff_sym(role) * raw"\, "
        sep = isempty(parts) ? (negative ? "-" : "") : (negative ? raw" - " : raw" + ")
        push!(parts, sep * csym * fsym)
        push!(sig, (role, negative ? -1 : 1, isunit))
    end
    return join(parts, ""), sig, used_symbol
end

_sig_string(sig) = join(["$(r):$(s):$(u)" for (r, s, u) in sort(sig)], ",")

const _COEFF_NOTE = "a_• are technology-specific conversion/efficiency coefficients (values differ by instance)"

# Build balance entries (key, formulation, element-kind) for a vertex.
function balance_entries(v::Node, asset)
    out = Tuple{String,ConstraintFormulation,String}[]
    for bid in balance_ids(v)
        bd = balance_data(v, bid)
        op = _balance_op(bd.sense)
        netflow = raw"\sum_{e \in \mathcal{I}(n)} f_{e,t} - \sum_{e \in \mathcal{O}(n)} f_{e,t}"
        syms = [raw"\mathcal{I}(n)" => "edges flowing into n", raw"\mathcal{O}(n)" => "edges flowing out of n"]
        if bid == :demand
            body = "$(netflow) + \\sum_{s} \\nu_{n,s,t} + \\sum_{s} q_{n,s,t} $(op) d_{n,t}"
            conds = ["ν, q present only where the node has non-served demand / supply segments"]
        else
            body = "$(netflow) $(op) 0"
            conds = String[]
        end
        f = ConstraintFormulation(name = "Node balance — :$(bid)", body = [body],
            indexsets = [raw"n \in \mathcal{N},\ t \in \mathcal{T}"], symbols = syms, conditions = conds)
        push!(out, ("BAL|Node|$(bid)|$(bd.sense)", f, "Node"))
    end
    return out
end

function balance_entries(v::AbstractStorage, asset)
    out = Tuple{String,ConstraintFormulation,String}[]
    kind = string(nameof(typeof(v)))
    for bid in balance_ids(v)
        bd = balance_data(v, bid)
        op = _balance_op(bd.sense)
        rhs, sig, used_symbol = _reconstruct_flow_terms(v, asset, bd)
        isempty(rhs) && (rhs = "0")
        conds = String[]
        if bid == :storage
            body = raw"s_{g,t} - (1 - \lambda_{g})\, s_{g,t-1} " * op * " " * rhs
            v isa LongDurationStorage && push!(conds, "across subperiod boundaries the previous level is adjusted by the storage-change variable Δ")
            name = "Storage balance (law of motion)"
        else
            body = rhs * " " * op * " 0"
            name = "Storage balance — :$(bid)"
        end
        used_symbol && push!(conds, _COEFF_NOTE)
        f = ConstraintFormulation(name = name, body = [body], indexsets = [raw"t \in \mathcal{T}"], conditions = conds)
        push!(out, ("BAL|$(kind)|$(bid)|$(_sig_string(sig))|$(bd.sense)", f, kind))
    end
    return out
end

function balance_entries(v::AbstractVertex, asset)   # Transformation and others
    out = Tuple{String,ConstraintFormulation,String}[]
    kind = asset === nothing ? string(nameof(typeof(v))) : string(nameof(typeof(asset)))
    for bid in balance_ids(v)
        bd = balance_data(v, bid)
        op = _balance_op(bd.sense)
        lhs, sig, used_symbol = _reconstruct_flow_terms(v, asset, bd)
        isempty(lhs) && continue
        body = lhs * " " * op * " 0"
        conds = used_symbol ? [_COEFF_NOTE] : String[]
        f = ConstraintFormulation(name = "$(kind) balance — :$(bid)", body = [body], indexsets = [raw"t \in \mathcal{T}"], conditions = conds)
        push!(out, ("BAL|$(kind)|$(bid)|$(_sig_string(sig))|$(bd.sense)", f, kind))
    end
    return out
end

# -----------------------------------------------------------------------------
# Collection: walk the system, grouping constraints
# -----------------------------------------------------------------------------

"""
    collect_active_formulations(system::System) -> OrderedDict

Walk a built `System` and collect every active constraint into an ordered
registry keyed so that each distinct constraint is recorded once with the list
of element ids it applies to. Balance constraints are reconstructed from
`balance_data` (so their signed flow terms are shown) and grouped by structure.

!!! note
    Call this *after* the model has been built: some constraints (e.g.
    `LongDurationStorageChangeConstraint`, `AgeBasedRetirementConstraint`) are
    pushed onto an element's `constraints` field during model generation.
"""
function collect_active_formulations(system::System, period::Int=period_index(system))
    reg = OrderedCollections.OrderedDict{String,ActiveConstraintEntry}()

    function add_entry!(key, f, element_kind, category, element_id)
        if !haskey(reg, key)
            reg[key] = ActiveConstraintEntry(formulation = f, element_kind = element_kind, category = category)
        end
        entry = reg[key]
        (element_id in entry.elements) || push!(entry.elements, element_id)
        (period in entry.periods) || push!(entry.periods, period)
    end

    function record!(y, asset)
        kind = string(nameof(typeof(y)))
        for c in all_constraints(y)
            if c isa BalanceConstraint
                for (key, f, klabel) in balance_entries(y, asset)
                    add_entry!(key, f, klabel, 2, id(y))
                end
            else
                add_entry!("C|$(typeof(c))|$(kind)", formulation(c, y), kind, _constraint_category(c), id(y))
            end
        end
        for f in structural_formulations(y)
            add_entry!("S|$(f.name)|$(kind)", f, kind, _STRUCTURAL_CATEGORY, id(y))
        end
    end

    for n in system.locations
        n isa AbstractVertex && record!(n, nothing)
    end
    for a in system.assets
        for f in fieldnames(typeof(a))
            y = getfield(a, f)
            (y isa AbstractEdge || y isa AbstractVertex) && record!(y, a)
        end
    end

    for (f, element_ids) in system_structural_formulations(system)
        for eid in element_ids
            add_entry!("SYS|$(f.name)", f, "system-wide", _STRUCTURAL_CATEGORY, eid)
        end
    end

    return reg
end

function _merge_registry!(into, from)
    for (key, entry) in from
        if haskey(into, key)
            target = into[key]
            for e in entry.elements
                (e in target.elements) || push!(target.elements, e)
            end
            for p in entry.periods
                (p in target.periods) || push!(target.periods, p)
            end
        else
            into[key] = entry
        end
    end
    return into
end

function collect_active_formulations(case::Case)
    reg = OrderedCollections.OrderedDict{String,ActiveConstraintEntry}()
    for system in get_periods(case)
        _merge_registry!(reg, collect_active_formulations(system, period_index(system)))
    end
    return reg
end

function collect_active_variables(system::System)
    seen = Set{String}()
    vars = VariableFormulation[]
    function record!(y)
        for v in formulation_variables(y)
            if !(v.symbol in seen)
                push!(seen, v.symbol)
                push!(vars, v)
            end
        end
    end
    for n in system.locations
        n isa MacroObject && record!(n)
    end
    for a in system.assets
        for f in fieldnames(typeof(a))
            y = getfield(a, f)
            y isa MacroObject && record!(y)
        end
    end
    return vars
end

collect_active_variables(case::Case) = collect_active_variables(first(get_periods(case)))

# -----------------------------------------------------------------------------
# Authored sections: sets and objective
# -----------------------------------------------------------------------------

function formulation_sets(system::System)
    lines = Pair{String,String}[]
    n_assets = length(system.assets); n_edges = length(get_edges(system))
    n_storages = length(get_storages(system)); n_nodes = length(get_nodes(system))
    td = isempty(system.time_data) ? nothing : first(values(system.time_data))
    n_t = td === nothing ? 0 : length(td.time_interval)
    n_w = td === nothing ? 0 : length(td.subperiod_indices)
    push!(lines, raw"\mathcal{S}" => "set of planning periods")
    push!(lines, raw"\mathcal{T}" => "set of modeled time steps ($(n_t) per representative horizon)")
    push!(lines, raw"\mathcal{W}" => "set of subperiods ($(n_w))")
    push!(lines, raw"\mathcal{A}" => "set of assets ($(n_assets))")
    push!(lines, raw"\mathcal{E}" => "set of edges ($(n_edges))")
    push!(lines, raw"\mathcal{G}" => "set of storages ($(n_storages))")
    push!(lines, raw"\mathcal{N}" => "set of nodes ($(n_nodes))")
    push!(lines, raw"\mathcal{I}(n),\ \mathcal{O}(n)" => "edges flowing into / out of node n")
    return lines
end

formulation_sets(case::Case) = formulation_sets(first(get_periods(case)))

function formulation_objective(system::System)
    has_uc = any(y -> y isa EdgeWithUC, get_edges(system))
    body = raw"\min \quad \sum_{s \in \mathcal{S}} \mathrm{df}_{s} \big( \mathrm{InvCost}_{s} + \mathrm{FOM}_{s} \big) + \sum_{s \in \mathcal{S}} \mathrm{df}_{s}\, \mathrm{opex}_{s}\, \mathrm{VarCost}_{s}"
    var = raw"\mathrm{VarCost}_{s} = \sum_{e, t} \omega_{t}\, c^{\mathrm{vom}}_{e}\, f_{e,t}"
    has_uc && (var *= raw" + \sum_{e \in \mathcal{E}^{UC}, t} c^{\mathrm{su}}_{e}\, K_{e}\, u^{+}_{e,t}")
    var *= raw" + (\text{policy slack penalties, if any})"
    symbols = Pair{String,String}[
        raw"\mathrm{df}_{s}" => "present-value discount factor for period s",
        raw"\mathrm{opex}_{s}" => "present-value annuity factor for operational costs in period s",
        raw"\mathrm{InvCost}_{s}" => raw"\sum_{y} c^{\mathrm{inv}}_{y}\, x^{\mathrm{new}}_{y}\ \text{(expandable edges/storages)}",
        raw"\mathrm{FOM}_{s}" => raw"\sum_{y} c^{\mathrm{fom}}_{y}\, x_{y}",
        raw"c^{\mathrm{vom}}_{e},\ c^{\mathrm{su}}_{e}" => "variable O&M and startup cost of edge e",
    ]
    return ConstraintFormulation(name = "Objective", body = [body, var], symbols = symbols)
end

formulation_objective(case::Case) = formulation_objective(first(get_periods(case)))

# -----------------------------------------------------------------------------
# Rendering helpers
# -----------------------------------------------------------------------------

_tex_escape(s::AbstractString) = replace(string(s),
    "_" => "\\_\\allowbreak ", "&" => "\\&", "%" => "\\%", "#" => "\\#",
    ">" => raw"$>$", "<" => raw"$<$",
    "≤" => raw"$\leq$", "≥" => raw"$\geq$", "×" => raw"$\times$",
    "↑" => raw"$\uparrow$", "↓" => raw"$\downarrow$", "→" => raw"$\rightarrow$",
    "Δ" => raw"$\Delta$", "ν" => raw"$\nu$", "ξ" => raw"$\xi$", "σ" => raw"$\sigma$",
    "₂" => raw"\textsubscript{2}")

_looks_like_math(s::AbstractString) = occursin('\\', s)
_md_desc(s::AbstractString) = _looks_like_math(s) ? "``$(s)``" : s
_tex_desc(s::AbstractString) = _looks_like_math(s) ? "\$" * s * "\$" : _tex_escape(s)

function _sorted_entries(reg)
    entries = collect(values(reg))
    sort!(entries, by = e -> (e.category, e.formulation.name))
    return entries
end

function _applies_to_md(entry::ActiveConstraintEntry; max_show::Int = 12)
    ids = string.(entry.elements)
    shown = length(ids) <= max_show ? ids : ids[1:max_show]
    s = join(["`$(x)`" for x in shown], ", ")
    length(ids) > max_show && (s *= ", … (+$(length(ids) - max_show) more)")
    return s
end

function _applies_to_tex(entry::ActiveConstraintEntry; max_show::Int = 12)
    ids = string.(entry.elements)
    shown = length(ids) <= max_show ? ids : ids[1:max_show]
    s = join([_tex_escape(x) for x in shown], ", ")
    length(ids) > max_show && (s *= ", \\ldots\\ (+$(length(ids) - max_show) more)")
    return s
end

# -----------------------------------------------------------------------------
# Markdown writer
# -----------------------------------------------------------------------------

function write_formulation_md(io::IO, sets, vars, params, objective, reg; multi_period::Bool=false)
    println(io, "# Model formulation\n")
    println(io, "_Auto-generated from the built model. Each constraint is shown once in generic form with the elements it applies to._\n")

    println(io, "## Sets and indices\n")
    for (sym, desc) in sets
        println(io, "- ``$(sym)``: $(desc)")
    end
    println(io)

    println(io, "## Decision variables\n")
    for v in vars
        dom = isempty(v.domain) ? "" : " (``$(v.domain)``)"
        idx = isempty(v.indexset) ? "" : ", for ``$(v.indexset)``"
        println(io, "- ``$(v.symbol)``$(dom): $(v.name)$(idx)")
    end
    println(io)

    println(io, "## Parameters\n")
    for (sym, desc) in params
        println(io, "- ``$(sym)``: $(desc)")
    end
    println(io)

    println(io, "## Objective\n")
    _md_constraint_block(io, objective; with_applies=false)

    cat_titles = Dict(1 => "Planning constraints", 2 => "Operation constraints", 3 => "Policy constraints", 4 => "Structural constraints (added automatically)")
    last_cat = 0
    println(io, "## Constraints\n")
    for entry in _sorted_entries(reg)
        if entry.category != last_cat
            println(io, "### $(get(cat_titles, entry.category, "Constraints"))\n")
            last_cat = entry.category
        end
        _md_constraint_block(io, entry.formulation; with_applies=true, entry=entry, multi_period=multi_period)
    end
    return nothing
end

function _md_constraint_block(io::IO, f::ConstraintFormulation; with_applies::Bool, entry=nothing, multi_period::Bool=false)
    println(io, "#### $(f.name)\n")
    if !isempty(f.body)
        println(io, "```math")
        println(io, "\\begin{aligned}")
        println(io, join(f.body, " \\\\\n"))
        println(io, "\\end{aligned}")
        println(io, "```\n")
    end
    if !isempty(f.indexsets)
        println(io, "for ", join(["``$(s)``" for s in f.indexsets], ", "), ".\n")
    end
    for (sym, desc) in f.symbols
        println(io, "- ``$(sym)``: $(_md_desc(desc))")
    end
    isempty(f.symbols) || println(io)
    for c in f.conditions
        println(io, "- _$(c)_")
    end
    isempty(f.conditions) || println(io)
    if with_applies && entry !== nothing
        println(io, "**Applies to** ($(entry.element_kind)): ", _applies_to_md(entry))
        multi_period && println(io, "  \n**Active in periods**: ", join(sort(entry.periods), ", "))
        println(io)
    end
    return nothing
end

# -----------------------------------------------------------------------------
# LaTeX writer
# -----------------------------------------------------------------------------

function write_formulation_tex(io::IO, sets, vars, params, objective, reg; multi_period::Bool=false)
    println(io, raw"\documentclass{article}")
    println(io, raw"\usepackage{amsmath,amssymb}")
    println(io, raw"\usepackage{graphicx}")
    println(io, raw"\usepackage[margin=1in]{geometry}")
    println(io, raw"\setlength{\parindent}{0pt}")
    println(io, raw"\sloppy")
    println(io, raw"\title{Model formulation}")
    println(io, raw"\date{}")
    println(io, raw"\begin{document}")
    println(io, raw"\maketitle")
    println(io, raw"\emph{Auto-generated from the built model. Each constraint is shown once in generic form with the elements it applies to.}")

    println(io, raw"\section*{Sets and indices}")
    println(io, raw"\begin{itemize}")
    for (sym, desc) in sets
        println(io, raw"\item $", sym, raw"$: ", _tex_escape(desc))
    end
    println(io, raw"\end{itemize}")

    println(io, raw"\section*{Decision variables}")
    println(io, raw"\begin{itemize}")
    for v in vars
        dom = isempty(v.domain) ? "" : " (\$" * v.domain * "\$)"
        idx = isempty(v.indexset) ? "" : ", for \$" * v.indexset * "\$"
        println(io, raw"\item $", v.symbol, raw"$", dom, ": ", _tex_escape(v.name), idx)
    end
    println(io, raw"\end{itemize}")

    println(io, raw"\section*{Parameters}")
    println(io, raw"\begin{itemize}")
    for (sym, desc) in params
        println(io, raw"\item $", sym, raw"$: ", _tex_desc(desc))
    end
    println(io, raw"\end{itemize}")

    println(io, raw"\section*{Objective}")
    _tex_constraint_block(io, objective; with_applies=false)

    cat_titles = Dict(1 => "Planning constraints", 2 => "Operation constraints", 3 => "Policy constraints", 4 => "Structural constraints (added automatically)")
    last_cat = 0
    println(io, raw"\section*{Constraints}")
    for entry in _sorted_entries(reg)
        if entry.category != last_cat
            println(io, raw"\subsection*{", get(cat_titles, entry.category, "Constraints"), "}")
            last_cat = entry.category
        end
        _tex_constraint_block(io, entry.formulation; with_applies=true, entry=entry, multi_period=multi_period)
    end
    println(io, raw"\end{document}")
    return nothing
end

function _tex_constraint_block(io::IO, f::ConstraintFormulation; with_applies::Bool, entry=nothing, multi_period::Bool=false)
    println(io, raw"\subsubsection*{", _tex_escape(f.name), "}")
    for eq in f.body
        println(io, raw"\begin{center}")
        println(io, raw"\resizebox{\ifdim\width>\linewidth \linewidth\else \width\fi}{!}{$\displaystyle ", eq, raw"$}")
        println(io, raw"\end{center}")
    end
    if !isempty(f.indexsets)
        println(io, "for ", join(["\$" * s * "\$" for s in f.indexsets], ", "), ".")
    end
    if !isempty(f.symbols)
        println(io, raw"\begin{itemize}")
        for (sym, desc) in f.symbols
            println(io, raw"\item $", sym, raw"$: ", _tex_desc(desc))
        end
        println(io, raw"\end{itemize}")
    end
    for c in f.conditions
        println(io, raw"\emph{", _tex_escape(c), "}\\\\")
    end
    if with_applies && entry !== nothing
        println(io, raw"\textbf{Applies to} (", _tex_escape(entry.element_kind), "): ", _applies_to_tex(entry), raw"\\")
        multi_period && println(io, raw"\textbf{Active in periods}: ", join(sort(entry.periods), ", "), raw"\\")
    end
    println(io)
    return nothing
end

# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

"""
    write_formulation(case::Case; path="formulation", formats=(:tex, :md))
    write_formulation(system::System; path="formulation", formats=(:tex, :md))

Write a compact, human-readable mathematical formulation of the model that was
built for `case`/`system`, reflecting only the constraints and variables that
are actually active. Produces a `.tex` file (compilable to PDF) and/or a `.md`
file (renderable in the Documenter docs via KaTeX), selected by `formats`.

`path` is the output path *without* extension; the relevant suffix is appended.

!!! note
    Call this *after* `generate_model` so that constraints added during model
    generation (e.g. long-duration storage and age-based retirement constraints)
    are included.
"""
function write_formulation(case::Case; path::AbstractString="formulation", formats=(:tex, :md))
    reg = collect_active_formulations(case)
    vars = collect_active_variables(case)
    sets = formulation_sets(case)
    params = formulation_parameters()
    obj = formulation_objective(case)
    multi_period = number_of_periods(case) > 1
    return _emit_formulation(path, formats, sets, vars, params, obj, reg, multi_period)
end

function write_formulation(system::System; path::AbstractString="formulation", formats=(:tex, :md))
    reg = collect_active_formulations(system)
    vars = collect_active_variables(system)
    sets = formulation_sets(system)
    params = formulation_parameters()
    obj = formulation_objective(system)
    return _emit_formulation(path, formats, sets, vars, params, obj, reg, false)
end

function _emit_formulation(path, formats, sets, vars, params, obj, reg, multi_period)
    written = String[]
    dir = dirname(path)
    isempty(dir) || mkpath(dir)
    if :md in formats
        md_path = endswith(path, ".md") ? path : path * ".md"
        open(md_path, "w") do io
            write_formulation_md(io, sets, vars, params, obj, reg; multi_period=multi_period)
        end
        @info "Wrote model formulation to $md_path"
        push!(written, md_path)
    end
    if :tex in formats
        tex_path = endswith(path, ".tex") ? path : path * ".tex"
        open(tex_path, "w") do io
            write_formulation_tex(io, sets, vars, params, obj, reg; multi_period=multi_period)
        end
        @info "Wrote model formulation to $tex_path"
        push!(written, tex_path)
    end
    return written
end
