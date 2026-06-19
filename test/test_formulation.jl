using Test
using HiGHS
import MacroEnergy:
    load_case,
    generate_model,
    create_optimizer,
    solution_algorithm,
    write_formulation,
    collect_active_formulations,
    formulation,
    AbstractTypeConstraint

const formulation_test_path = joinpath(@__DIR__, "test_inputs")

function test_formulation_outputs()
    case = load_case(formulation_test_path)
    optimizer = create_optimizer(HiGHS.Optimizer)
    alg = solution_algorithm(case)
    model = generate_model(case, optimizer, alg)

    outdir = mktempdir()
    base = joinpath(outdir, "formulation")
    written = write_formulation(case; path = base)

    md_path = base * ".md"
    tex_path = base * ".tex"

    # Both files are produced and non-empty
    @test Set(written) == Set([md_path, tex_path])
    @test isfile(md_path) && filesize(md_path) > 0
    @test isfile(tex_path) && filesize(tex_path) > 0

    md = read(md_path, String)
    tex = read(tex_path, String)

    # Expected document sections
    for section in ("Sets and indices", "Decision variables", "Parameters", "Objective", "Constraints")
        @test occursin(section, md)
    end
    @test occursin(raw"\documentclass{article}", tex)
    @test occursin(raw"\end{document}", tex)

    # Compact-symbol constraint math
    @test occursin(raw"f_{e,t} \leq \alpha_{e,t}\, x_{e}", md)
    @test occursin("Applies to", md)

    # Balances are reconstructed with their signed flow terms (issue 1b)
    @test occursin("Node balance", md)
    @test occursin(raw"\sum_{e \in \mathcal{I}(n)} f_{e,t} - \sum_{e \in \mathcal{O}(n)} f_{e,t}", md)
    # Storage law of motion shows the actual signed flows
    @test occursin(raw"s_{g,t} - (1 - \lambda_{g})\, s_{g,t-1}", md)
    # Reconstructed coefficients are symbolic (a_role), not numeric
    @test occursin(raw"a_{\text{", md)
    @test !occursin(r"\d\.\d+\\, f_", md)   # no "1.234\, f_..." numeric coefficients

    # Structural constraints (added with @constraint, not via add_model_constraint!)
    @test occursin("Structural constraints", md)
    @test occursin(raw"x_{e} = x^{\mathrm{new}}_{e} - x^{\mathrm{ret}}_{e} + x^{\mathrm{ex}}_{e}", md)
    @test occursin(raw"u_{e,t} - u_{e,t-1} = u^{+}_{e,t} - u^{-}_{e,t}", md)

    # Every active constraint must resolve to an authored formulation (no fallback)
    reg = collect_active_formulations(case)
    @test !isempty(reg)
    fallback_hits = [
        (key, entry) for (key, entry) in reg
        if any(c -> occursin("not yet authored", c), entry.formulation.conditions)
    ]
    @test isempty(fallback_hits)

    return nothing
end

@testset "Model formulation export" begin
    test_formulation_outputs()
end
