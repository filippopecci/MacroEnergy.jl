function add_kirchoff_voltage_law_constraint!(system::System, model::Model)
    
    all_power_lines = [a.elec_edge for a in system.assets if a isa PowerLine]
    kvl_power_lines = [e for e in all_power_lines if existing_capacity(e)>0 && !ismissing(e.kvl_coefficient)]

    edge_list = [(findfirst(n.id==e.start_vertex.id for n in system.locations),
                 findfirst(n.id==e.end_vertex.id for n in system.locations)) for e in kvl_power_lines]
    if !isempty(edge_list)
        cycle_incidence = get_cycle_basis_incidence(edge_list)
        @info " -- Adding Kirchoff Voltage Law constraint"
        if isempty(cycle_basis_incidence)
            @warn " -- The subgraph of power lines subject to KVL does not have cycles, no constraint will added"
        else
            @constraint(model, 
                [c in axes(cycle_incidence,1)],
                sum(cycle_basis_incidence[c,i] * kvl_coefficient(kvl_power_lines[i]) * flow(kvl_power_lines[i]) for i in eachindex(kvl_power_lines)) == 0
            )
        end
    end

    return nothing
end

function get_cycle_basis_incidence(elist::Vector{Tuple{Int, Int}})

    cycle_basis = Graphs.cycle_basis(Graphs.SimpleGraph(Graphs.SimpleEdge.(elist)))
    # Construct cycle basis incidence matrix
    num_cycles = length(cycle_basis)
    cycle_incidence = zeros(Int, num_cycles, length(elist))

    for (ci, cycle) in enumerate(cycle_basis)
        # Convert cycle to a list of edges (as pairs of node indices)
        cycle_edges = [(cycle[i], cycle[mod1(i+1, length(cycle))]) for i in eachindex(cycle)]
        # For each edge in the cycle, find its index in elist (or reversed)
        for (u, v) in cycle_edges
            idx = findfirst(e -> (e == (u, v)) || (e == (v, u)), elist)
            if idx !== nothing
                # Assign +1 if edge orientation matches, -1 otherwise
                cycle_incidence[ci, idx] = elist[idx] == (u, v) ? 1 : -1
            end
        end
    end

    return cycle_incidence
end