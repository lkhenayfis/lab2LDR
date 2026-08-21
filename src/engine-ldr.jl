
"""
    lab2LDR(path::String)

Le a configuracao de estudo do lab em `path` e monta um objeto LDRModel correspondente
"""
function lab2LDR(path::String, optimizer)::LDRModel
    study = SDDPlab.read_study(path)
    files = study.inputs.files

    system = SDDPlab.get_system(files)
    scenarios = SDDPlab.get_scenarios(files)

    num_stages = get_number_of_stages(get_graph(scenarios))

    m = LDRModel(optimizer)

    add_system_elements!(m, system, num_stages)
    add_hydro_balance!(m, get_hydros(system))
end

# STAGE ATTRIBUTE ----------------------------------------------------------------------------------

function build_cartesian_indexes(a::Array)
    dim = ndims(a)
    idx = [axes(a, i) for i in 1:(dim - 1)]
    return CartesianIndex.(idx...)
end

function set_stage_attribute!(m::LDRModel, v::Symbol, num_stages::Int)
    cart_indexes = build_cartesian_indexes(m[v])
    for s in 1:num_stages
        set_attribute.(m[v][cart_indexes, s], [Stage(s)])
    end 
end

function set_stage_attribute!(m::LDRModel, vs::Vector{Symbol}, num_stages::Int)
    for v in vs
        cart_indexes = build_cartesian_indexes(m[v])
        for s in 1:num_stages
            set_attribute.(m[v][cart_indexes, s], [Stage(s)])
        end
    end
end

# SYSTEM ELEMENT ADDERS ----------------------------------------------------------------------------

function add_system_elements!(m::LDRModel, ses::Buses, num_stages::Int)
    num_buses = length(ses)

    m[LOAD] = JuMP.@variable(m, [1:num_buses, 1:num_stages], base_name = String(LOAD))
    m[DEFICIT] = JuMP.@variable(m, [1:num_buses, 1:num_stages], base_name = String(DEFICIT))
    
    for n in 1:num_buses, s in 1:num_stages
        JuMP.set_lower_bound(m[DEFICIT][n, s], 0)
    end

    set_stage_attribute!(m, [LOAD, DEFICIT], num_stages)

    return nothing
end

function add_system_elements!(m::LDRModel, ses::Lines, num_stages::Int)
    num_lines = length(ses)

    m[DIRECT_EXCHANGE] = JuMP.@variable(
        m, [1:num_lines, 1:num_stages], base_name = String(DIRECT_EXCHANGE)
    )
    for n in 1:num_lines, s in 1:num_stages
        JuMP.set_lower_bound(m[DIRECT_EXCHANGE][n, s], 0)
        JuMP.set_upper_bound(m[DIRECT_EXCHANGE][n, s], ses.entities[n].capacity)
    end

    m[REVERSE_EXCHANGE] = JuMP.@variable(
        m, [1:num_lines, 1:num_stages], base_name = String(REVERSE_EXCHANGE)
    )
    for n in 1:num_lines, s in 1:num_stages
        JuMP.set_lower_bound(m[REVERSE_EXCHANGE][n, s], 0)
        JuMP.set_upper_bound(m[REVERSE_EXCHANGE][n, s], ses.entities[n].capacity)
    end

    return m[NET_EXCHANGE] = JuMP.@expression(
        m,
        [n = 1:num_lines, 1:num_stages],
        m[DIRECT_EXCHANGE][n, s] - m[REVERSE_EXCHANGE][n, s]
    )

    set_stage_attribute!(m, [DIRECT_EXCHANGE, REVERSE_EXCHANGE, NET_EXCHANGE], num_blocks)

    return nothing
end

function add_system_elements!(m::LDRModel, ses::Thermals, num_stages::Int)
    num_thermals = length(ses)
    m[THERMAL_GENERATION] = JuMP.@variable(
        m, [1:num_thermals, 1:num_stages], base_name = String(THERMAL_GENERATION)
    )
    for n in 1:num_thermals, s in 1:num_stages
        JuMP.set_lower_bound(m[THERMAL_GENERATION][n, s], ses.entities[n].min_generation)
        JuMP.set_upper_bound(m[THERMAL_GENERATION][n, s], ses.entities[n].max_generation)
    end

    m[THERMAL_GENERATION_COST] = JuMP.@expression(
        m,
        [n = 1:num_thermals, s = 1:num_stages],
        ses.entities[n].cost * m[THERMAL_GENERATION][n, s]
    )

    set_stage_attribute!(m, THERMAL_GENERATION, num_stages)

    return nothing
end

function add_system_elements!(m::LDRModel, ses::Hydros, num_stages::Int)
    num_hydros = length(ses)

    m[STORED_VOLUME] = JuMP.@variable(
        m, [n = 1:num_hydros, s = 1:num_stages], base_name = String(STORED_VOLUME),
    )

    for n in 1:num_hydros, s in 1:num_stages
        JuMP.set_lower_bound(m[STORED_VOLUME][n, s], ses.entities[n].min_storage)
        JuMP.set_upper_bound(m[STORED_VOLUME][n, s], ses.entities[n].max_storage)
    end

    # dimensao 2 gambiarrada com comprimento 1 para manter consistencia com todo o resto em
    # set_stage_attribute!
    m[INFLOW] = JuMP.@variable(m, [1:num_hydros, 1:num_stages], base_name = String(INFLOW))

    m[TURBINED_FLOW] = JuMP.@variable(
        m, [1:num_hydros, 1:num_stages], base_name = String(TURBINED_FLOW)
    )
    for n in 1:num_hydros, s in 1:num_stages
        JuMP.set_lower_bound(
            m[TURBINED_FLOW][n, s],
            ses.entities[n].min_generation / ses.entities[n].productivity,
        )
        JuMP.set_upper_bound(
            m[TURBINED_FLOW][n, s],
            ses.entities[n].max_generation / ses.entities[n].productivity,
        )
    end

    m[SPILLAGE] = JuMP.@variable(
        m, [1:num_hydros, 1:num_stages], base_name = String(SPILLAGE)
    )
    for n in 1:num_hydros, s in 1:num_stages
        JuMP.set_lower_bound(m[SPILLAGE][n, s], 0)
    end

    m[OUTFLOW] = JuMP.@expression(
        m, [n = 1:num_hydros, s = 1:num_stages], m[TURBINED_FLOW][n, s] + m[SPILLAGE][n, s]
    )

    m[HYDRO_GENERATION] = JuMP.@expression(
        m,
        [n = 1:num_hydros, s = 1:num_stages],
        ses.entities[n].productivity * m[TURBINED_FLOW][n, s]
    )

    m[HYDRO_MIN_GENERATION_SLACK] = JuMP.@variable(
        m, [1:num_hydros, 1:num_stages], base_name = String(HYDRO_MIN_GENERATION_SLACK)
    )
    for n in 1:num_hydros, s in 1:num_stages
        JuMP.set_lower_bound(m[HYDRO_MIN_GENERATION_SLACK][n, s], 0)
    end

    JuMP.@constraint(
        m,
        [n = 1:num_hydros, s = 1:num_stages],
        m[HYDRO_GENERATION][n, s] + m[HYDRO_MIN_GENERATION_SLACK][n, s] >=
            ses.entities[n].min_generation
    )

    set_stage_attribute!(
        m,
        [STORED_VOLUME, INFLOW, TURBINED_FLOW, SPILLAGE, HYDRO_MIN_GENERATION_SLACK],
        num_stages
    )

    return nothing
end

function add_system_elements!(m::LDRModel, ses::SystemData, num_stages::Int)
    add_system_elements!(m, get_buses(ses), num_stages)
    add_system_elements!(m, get_lines(ses), num_stages)
    add_system_elements!(m, get_thermals(ses), num_stages)
    add_system_elements!(m, get_hydros(ses), num_stages)
    return nothing
end

# HYDRO BALANCE ------------------------------------------------------------------------------------

function add_hydro_balance!(m::LDRModel, hydros::Hydros)
    num_hydros = length(hydros)

    m[HYDRO_BALANCE] = JuMP.@constraint(
        m,
        [n = 1:num_hydros, s = 1:num_stages],
        m[STORED_VOLUME][n, s] ==
            (s > 1 ? m[STORED_VOLUME][n, s-1] : hydros.entities[n].initial_storage) - m[OUTFLOW][n] +
        m[INFLOW][n, s] +
        sum(
            m[OUTFLOW][j, s] for j in 1:num_hydros if
            downstream(hydros.entities[j].id, hydros) == hydros.entities[n]
        )
    )
    return nothing
end