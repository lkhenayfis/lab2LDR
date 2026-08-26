
"""
    lab2LDR(path::String)

Le a configuracao de estudo do lab em `path` e monta um objeto LDRModel correspondente
"""
function lab2LDR(path::String, optimizer)::LDRModel
    study = SDDPlab.read_study(path)
    files = study.inputs.files

    system = SDDPlab.get_system(files)
    scenarios = SDDPlab.get_scenarios(files)
    stochastic = SDDPlab.get_stochastic_process(scenarios.inflow)

    num_stages = get_number_of_stages(get_graph(scenarios))

    m = LDRModel(optimizer)

    add_system_elements!(m, system, num_stages)
    add_hydro_balance!(m, get_hydros(system))
    add_load_balance!(m, files, num_stages)
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

    #m[INFLOW] = JuMP.@variable(m, [1:num_hydros, 1:num_stages], base_name = String(INFLOW))

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
        [STORED_VOLUME, TURBINED_FLOW, SPILLAGE, HYDRO_MIN_GENERATION_SLACK],
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

# BALANCE CONSTRAINTS ------------------------------------------------------------------------------

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

function add_load_balance!(m::LDRModel, files::Vector{InputModule}, num_stages::Int)
    system = get_system(files)
    hydros_entities = get_hydros_entities(system)
    thermals_entities = get_thermals_entities(system)
    lines_entities = get_lines_entities(system)
    scenarios = get_scenarios(files)
    bus_ids = get_ids(get_buses(system))

    num_buses = length(bus_ids)
    num_lines = length(lines_entities)
    num_hydros = length(hydros_entities)
    num_thermals = length(thermals_entities)

    m[LOAD_BALANCE] = JuMP.@constraint(
        m,
        [n = 1:num_buses, s = 1:num_stages],
        sum(
            m[HYDRO_GENERATION][j, s] for
            j in 1:num_hydros if hydros_entities[j].bus_id == bus_ids[n]
        ) +
        sum(
            m[THERMAL_GENERATION][j, s] for
            j in 1:num_thermals if thermals_entities[j].bus_id == bus_ids[n]
        ) +
        sum(
            m[DIRECT_EXCHANGE][j, s] - m[REVERSE_EXCHANGE][j, s] for
            j in 1:num_lines if lines_entities[j].target_bus_id == bus_ids[n]
        ) +
        sum(
            m[REVERSE_EXCHANGE][j, s] - m[DIRECT_EXCHANGE][j, s] for
            j in 1:num_lines if lines_entities[j].source_bus_id == bus_ids[n]
        ) +
        m[DEFICIT][n, s] == get_load(bus_ids[n], s, scenarios)
    )
    return nothing
end

# UNCERTAINTY --------------------------------------------------------------------------------------

function get_models_per_stage(stch::Naive)::Vector{Vector{T} where T <: Distribution}
    num_hydros, num_stages = size(stch)
    models = []
    for i in 1:num_stages
        stage_i = [stch.models[j].distributions[i] for j in 1:num_hydros]
        push!(models, stage_i)
    end
    
    return models
end

function get_copulas_per_stage(stch::Naive)::Vector{Copula}
    [stch.copulas[j] for j in range(1, size(stch, 2))]
end

function validate_univar_gaussian(models_per_stage)
    valid = true
    for stage_i in models_per_stage
        for dist in stage_i
            valid = valid & (isa(dist, Normal) | isa(dist, LogNormal))
        end
    end
    
    if !valid
        error("todas as distribuicoes univariadas devem ser Normais ou LogNormais")
    else
        return true
    end
end

function validate_copulas_kind(copulas_per_stage)
    valid = true
    for copula in copulas_per_stage
        valid = valid & (isa(copula, IndependentCopula) | isa(copula, GaussianCopula))
    end

    if !valid
        error("todas as copulas devem ser Independentes ou Gaussianas")
    else
        return true
    end
end

function validate_stch_for_ldr(stch::Naive)
    models = get_models_per_stage(stch)
    copulas = get_copulas_per_stage(stch)
    aux = validate_univar_gaussian(models) && validate_copulas_kind(copulas)
    return models, copulas
end

function copula2Σ(c::IndependentCopula)
    diagm(repeat([0], length(c)))
end

function copula2Σ(c::GaussianCopula)
    Σ = c.Σ
    Σ[diagind(Σ)] .= 0.0
    return Σ
end

function build_stage_ldruncertainty_1d(
    stage_models::Vector{T} where T <: Distribution,
    stage_copula::Copula;
    truncate_mode = [:value, :quantile],
    truncate_at = [0, .975]
    )::UnivariateDistribution

    mod = stage_models[1]

    truncate_at[1] = (truncate_mode[1] == :quantile) ? quantile(mod, truncate_at[1]) : truncate_at[1]
    truncate_at[2] = (truncate_mode[2] == :quantile) ? quantile(mod, truncate_at[2]) : truncate_at[2]
    mod = truncated(mod, truncate_at...)

    return mod
end

function build_stage_ldruncertainty_Nd(
    stage_models::Vector{T} where T <: Distribution,
    stage_copula::IndependentCopula
    )::MultivariateDistribution

    return Distributions.product_distribution(stage_models...)    
end

function build_stage_ldruncertainty_Nd(
    stage_models::Vector{T} where T <: Distribution,
    stage_copula::Copula
    )::MultivariateDistribution

    means = [m.μ for m in stage_models]
    stds = [m.σ for m in stage_models]

    Σ = diagm(stds) + copula2Σ(stage_copula)

    if isa(stage_models[1], Normal)
        return MvNormal(means, Σ)
    elseif isa(stage_models[1], LogNormal)
        return MvLogNormal(means, Σ)
    end
end

function build_stage_ldruncertainty(
    stage_models::Vector{T} where T <: Distribution,
    stage_copula::Copula
    )::UnivariateDistribution

    if length(stage_models) == 1
        build_stage_ldruncertainty_1d(stage_models, stage_copula)
    else
        build_stage_ldruncertainty_Nd(stage_models, stage_copula)
    end
end

function stch2ldruncertainty(stch::Naive)::Vector{Distribution}
    models, copulas = validate_stch_for_ldr(stch)

    mvdists = []
    for (m,c) in zip(models, copulas)
        push!(mvdists, build_stage_ldruncertainty(m, c))
    end

    return mvdists
end

function determinisitic_stage1(dist::Distribution)::Vector{Float64}
    isa(dist, UnivariateDistribution) ? [mean(dist)] : mean(dist)
end

function add_inflow!(m::LDRModel, hydros::Hydros, stch::Naive)
    num_hydros, num_stages = size(stch)

    ldr_dists = stch2ldruncertainty(stch)
    stage_1 = determinisitic_stage1(ldr_dists[1])

    m[:INFLOW_1] = JuMP.@variable(m, [1:num_hydros], base_name = "INFLOW_1")
    for n in 1:num_hydros
        JuMP.fix(m[:INFLOW_1][n], stage_1[n])
    end

    for s in 2:num_stages
        add_uncertain_inflow!(m, s, ldr_dists[s])
    end
end

function add_uncertain_inflow!(m::LDRModel, stage::Int, dist::UnivariateDistribution)
    varname = Symbol("INFLOW_" * string(stage))
    JuMP.@variable(
        m,
        base_name = String(varname),
        set = Uncertainty(; distribution = dist)
    )
    return nothing
end

function add_uncertain_inflow!(m::LDRModel, stage::Int, dist::ProductDistribution)
    size = length(dist)
    varname = Symbol("INFLOW_" * string(stage))
    JuMP.@variable(
        m,
        [1:size],
        base_name = String(varname),
        set = Uncertainty(; distribution = dist)
    )
    return nothing
end

function add_uncertain_inflow!(m::LDRModel, stage::Int, dist::MultivariateNormal; alpha = .95)
    size = length(dist)
    dist = LinearDecisionRules.ConfidenceMvNormal(dist.μ, dist.Σ, alpha)
    varname = Symbol("INFLOW_" * string(stage))
    JuMP.@variable(
        m,
        [1:size],
        base_name = String(varname),
        set = Uncertainty(; distribution = dist)
    )
    return nothing
end