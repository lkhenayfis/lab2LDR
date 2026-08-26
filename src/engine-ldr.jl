
"""
    build_model(path::String, optimizer; n_breakpoints = 0)::LDRModel

Le a configuracao de estudo do lab (formato SDDPlab) em `path` e monta o `LDRModel`
correspondente, pronto para ser otimizado.

Esta e a funcao de entrada da ponte entre o SDDPlab e o LinearDecisionRules. Ela cumpre
o mesmo papel de `Lab.build` no motor SDDP do laboratorio, porem montando o problema
multiestagio como um modelo de regras de decisao lineares (LDR) ao inves de um modelo SDDP.

Fluxo:
1. `read_study` carrega o `main.jsonc` e demais arquivos de entrada;
2. as variaveis de sistema (cargas, deficit, intercambios, termicas e hidros) sao criadas
   por `add_system_elements!`;
3. `add_inflow!` declara as variaveis de afluencia — deterministica no primeiro estagio e
   incerta (`Uncertainty`) nos demais;
4. `add_hydro_balance!` e `add_load_balance!` adicionam as restricoes de balanco hidrico e
   de atendimento a carga;
5. `add_objective!` define a funcao objetivo de minimizacao de custo.

# Argumentos
- `path`: diretorio do estudo contendo `main.jsonc`.
- `optimizer`: fabrica de otimizador JuMP (ex.: `HiGHS.Optimizer`).
- `n_breakpoints`: numero de pontos de quebra interiores para aproximar a regra de decisao
  de cada afluencia incerta por partes (piecewise linear). `0` mantem a regra puramente afim.

# Observacao
Atualmente apenas o caso com uma unica usina hidreletrica esta funcional; os demais caminhos
(multiplas hidros, linhas de intercambio) ainda nao estao operacionais.
"""
function build_model(path::String, optimizer; n_breakpoints = 0)::LDRModel
    cuts = EndCuts(path)
    study = read_study(path)
    do_truncation!(path, study)

    files = study.inputs.files

    system = get_system(files)
    scenarios = get_scenarios(files)
    stochastic = scenarios.inflow.stochastic_process

    num_stages = get_number_of_stages(get_graph(scenarios))

    m = LDRModel(optimizer)

    add_system_elements!(m, system, num_stages)
    add_inflow!(m, stochastic; n_breakpoints = n_breakpoints)
    add_hydro_balance!(m, get_hydros(system), stochastic)
    add_load_balance!(m, files, num_stages)
    add_endcuts!(m, cuts)
    add_objective!(m, system, stochastic)

    return m
end

# STAGE ATTRIBUTE ----------------------------------------------------------------------------------

"""
    build_cartesian_indexes(a::Array)

Retorna o vetor de `CartesianIndex` que percorre todas as dimensoes de `a` exceto a ultima.

E usada por `set_stage_attribute!` para iterar sobre os indices "de entidade" de um container
de variaveis `[entidade..., estagio]`, mantendo a dimensao de estagio livre para ser fixada
externamente.
"""
function build_cartesian_indexes(a::Array)
    dim = ndims(a)
    idx = [axes(a, i) for i in 1:(dim - 1)]
    return CartesianIndex.(idx...)
end

"""
    set_stage_attribute!(m::LDRModel, v::Symbol, num_stages::Int)

Marca todas as variaveis registradas em `m` sob a chave `v` com o atributo `Stage` do
LinearDecisionRules, associando cada coluna (ultima dimensao do container) ao seu estagio.

O `Stage` informa ao LDR em que estagio cada decisao e tomada, o que define quais variaveis
sao nao-antecipativas (estagio 1) e quais podem adaptar-se afimente a incerteza revelada.
"""
function set_stage_attribute!(m::LDRModel, v::Symbol, num_stages::Int)
    cart_indexes = build_cartesian_indexes(m[v])
    for s in 1:num_stages
        set_attribute.(m[v][cart_indexes, s], [Stage(s)])
    end
end

"""
    set_stage_attribute!(m::LDRModel, vs::Vector{Symbol}, num_stages::Int)

Aplica [`set_stage_attribute!`](@ref) a cada chave de variavel em `vs`, marcando o atributo
`Stage` para varios grupos de variaveis de uma vez.
"""
function set_stage_attribute!(m::LDRModel, vs::Vector{Symbol}, num_stages::Int)
    for v in vs
        cart_indexes = build_cartesian_indexes(m[v])
        for s in 1:num_stages
            set_attribute.(m[v][cart_indexes, s], [Stage(s)])
        end
    end
end

# SYSTEM ELEMENT ADDERS ----------------------------------------------------------------------------

"""
    add_system_elements!(m::LDRModel, ses::Buses, num_stages::Int)

Adiciona ao modelo as variaveis associadas aos barramentos: a carga (`LOAD`) e o deficit de
atendimento (`DEFICIT`, restrito a valores nao-negativos), indexadas por barra e estagio, e
registra o atributo de estagio de ambas.
"""
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

"""
    add_system_elements!(m::LDRModel, ses::Lines, num_stages::Int)

Adiciona as variaveis de intercambio direto e reverso (`DIRECT_EXCHANGE`, `REVERSE_EXCHANGE`),
limitadas pela capacidade de cada linha, e define a expressao de intercambio liquido
(`NET_EXCHANGE`).

!!! warning
    Caminho ainda nao funcional. O codigo apos o primeiro `return` (incluindo
    `set_stage_attribute!`) e inalcancavel e a variavel `num_blocks` nao esta definida neste
    escopo. Mantido como esqueleto ate que o suporte a linhas seja implementado.
"""
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

"""
    add_system_elements!(m::LDRModel, ses::Thermals, num_stages::Int)

Adiciona as variaveis de geracao termica (`THERMAL_GENERATION`), limitadas pelos despachos
minimo e maximo de cada usina, a expressao de custo de geracao (`THERMAL_GENERATION_COST`) e
o atributo de estagio.
"""
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

"""
    add_system_elements!(m::LDRModel, ses::Hydros, num_stages::Int)

Adiciona as variaveis hidraulicas: volume armazenado (`STORED_VOLUME`, entre os limites de
armazenamento), vazao turbinada (`TURBINED_FLOW`, limitada pela geracao min/max dividida pela
produtividade) e vertimento (`SPILLAGE`, nao-negativo). Define ainda as expressoes de vazao
defluente total (`OUTFLOW = TURBINED_FLOW + SPILLAGE`) e de geracao hidraulica
(`HYDRO_GENERATION = produtividade * TURBINED_FLOW`), e registra o atributo de estagio das
variaveis de decisao.
"""
function add_system_elements!(m::LDRModel, ses::Hydros, num_stages::Int)
    num_hydros = length(ses)

    m[STORED_VOLUME] = JuMP.@variable(
        m, [n = 1:num_hydros, s = 1:num_stages], base_name = String(STORED_VOLUME),
    )

    for n in 1:num_hydros, s in 1:num_stages
        JuMP.set_lower_bound(m[STORED_VOLUME][n, s], ses.entities[n].min_storage)
        JuMP.set_upper_bound(m[STORED_VOLUME][n, s], ses.entities[n].max_storage)
    end

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

    set_stage_attribute!(
        m,
        [STORED_VOLUME, TURBINED_FLOW, SPILLAGE],
        num_stages
    )

    return nothing
end

"""
    add_system_elements!(m::LDRModel, ses::SystemData, num_stages::Int)

Conveniencia que adiciona ao modelo, em ordem, todos os grupos de elementos do sistema
(barras, linhas, termicas e hidros) delegando aos metodos especializados de
[`add_system_elements!`](@ref).
"""
function add_system_elements!(m::LDRModel, ses::SystemData, num_stages::Int)
    add_system_elements!(m, get_buses(ses), num_stages)
    add_system_elements!(m, get_lines(ses), num_stages)
    add_system_elements!(m, get_thermals(ses), num_stages)
    add_system_elements!(m, get_hydros(ses), num_stages)
    return nothing
end

# BALANCE CONSTRAINTS ------------------------------------------------------------------------------

"""
    add_hydro_balance!(m::LDRModel, hydros::Hydros, stch::Naive)

Adiciona a restricao de balanco hidrico (`HYDRO_BALANCE`) para cada usina e estagio:

    volume[s] == volume_inicial[s-1] - defluencia + afluencia + defluencia_de_montante

onde o volume do estagio anterior e substituido pelo armazenamento inicial no primeiro
estagio, a afluencia e a variavel `INFLOW_s` correspondente ao estagio, e a soma de montante
percorre as usinas cuja jusante (`downstream`) e a usina corrente.

O ramo `num_hydros == 1` trata a afluencia como escalar (exigencia do LinearDecisionRules para
incerteza univariada); o ramo com multiplas usinas usa afluencias vetoriais e ainda nao esta
operacional.
"""
function add_hydro_balance!(m::LDRModel, hydros::Hydros, stch::Naive)
    num_hydros = length(hydros)
    num_stages = size(stch, 2)

    inflow_vars = [Symbol("INFLOW_" * string(s)) for s in 1:num_stages]

    # gambiarra para lidar com INFLOW_t tendo tipos diferentes
    # no caso de incerteza univariada, a declaracao da variavel de inflow **PRECISA** ser escalar
    # do contrario da erro do LinearDecisionRules
    # sendo assim, quando tem so uma hidro todos os INFLOW_t sao escalares
    # quando tem mais de uma elas sao vetores com comprimento num_stages
    if num_hydros == 1
        m[HYDRO_BALANCE] = JuMP.@constraint(
            m,
            [s = 1:num_stages],
            m[STORED_VOLUME][1, s] ==
                (s > 1 ? m[STORED_VOLUME][1, s-1] : hydros.entities[1].initial_storage) - m[OUTFLOW][1] +
            m[inflow_vars[s]] +
            sum(
                m[OUTFLOW][j, s] for j in 1:num_hydros if
                downstream(hydros.entities[j].id, hydros) == hydros.entities[1]
            )
        )
    else
        m[HYDRO_BALANCE] = JuMP.@constraint(
            m,
            [n = 1:num_hydros, s = 1:num_stages],
            m[STORED_VOLUME][n, s] ==
                (s > 1 ? m[STORED_VOLUME][n, s-1] : hydros.entities[n].initial_storage) - m[OUTFLOW][n] +
            m[inflow_vars[s]][n] +
            sum(
                m[OUTFLOW][j, s] for j in 1:num_hydros if
                downstream(hydros.entities[j].id, hydros) == hydros.entities[n]
            )
        )
    end
    return nothing
end

"""
    add_load_balance!(m::LDRModel, files::Vector{InputModule}, num_stages::Int)

Adiciona a restricao de atendimento a carga (`LOAD_BALANCE`) por barra e estagio: a soma da
geracao hidraulica e termica conectada a barra, do intercambio liquido que chega/sai e do
deficit deve igualar a carga daquela barra e estagio (`get_load`).
"""
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

"""
    get_models_per_stage_1d(stch::Naive)::Vector{<:Distribution}

Extrai, para o caso de uma unica usina, o vetor de distribuicoes marginais de afluencia
indexado por estagio a partir do processo estocastico `Naive`.
"""
function get_models_per_stage_1d(stch::Naive)::Vector{<:Distribution}
    num_stages = size(stch, 2)
    models = Vector{Distribution}()
    for i in 1:num_stages
        stage_i = stch.models[1].distributions[i]
        push!(models, stage_i)
    end
    return models
end

"""
    get_models_per_stage_Nd(stch::Naive)::Vector{Vector{<:Distribution}}

Extrai, para o caso multivariado, as distribuicoes marginais de afluencia organizadas como um
vetor por estagio, cada um contendo a marginal de cada usina.

!!! warning
    Caminho ainda nao funcional: a variavel `num_hydros` nao esta definida neste escopo.
"""
function get_models_per_stage_Nd(stch::Naive)::Vector{Vector{<:Distribution}}
    num_stages = size(stch, 2)
    models = Vector{Vector{Distribution}}()
    for i in 1:num_stages
        stage_i = [stch.models[j].distributions[i] for j in 1:num_hydros]
        push!(models, stage_i)
    end
    return models
end

"""
    get_models_per_stage(stch::Naive)

Despacha entre [`get_models_per_stage_1d`](@ref) e [`get_models_per_stage_Nd`](@ref) conforme o
numero de usinas (`size(stch, 1)`), retornando as distribuicoes marginais de afluencia por
estagio.
"""
function get_models_per_stage(stch::Naive)::Union{Vector{<:Distribution},Vector{Vector{<:Distribution}}}
    num_hydros = size(stch, 1)

    if num_hydros == 1
        models = get_models_per_stage_1d(stch)
    else
        models = get_models_per_stage_Nd(stch)
    end

    return models
end

"""
    get_copulas_per_stage(stch::Naive)::Vector{Copula}

Retorna o vetor de copulas do processo `Naive` indexado por estagio, usadas para descrever a
dependencia entre as afluencias das usinas em cada estagio.
"""
function get_copulas_per_stage(stch::Naive)::Vector{Copula}
    [stch.copulas[j] for j in range(1, size(stch, 2))]
end

"""
    validate_univar_gaussian(model_per_stage::Vector{<:Distribution})

Garante que todas as marginais (caso univariado) sejam `Normal` ou `LogNormal` — as unicas
familias suportadas pela reformulacao LDR desta ponte. Lanca erro caso contrario, retorna
`true` se validas.
"""
function validate_univar_gaussian(model_per_stage::Vector{<:Distribution})
    valid = true
    for dist in model_per_stage
        valid = valid & (isa(dist, Normal) | isa(dist, LogNormal))
    end

    if !valid
        error("todas as distribuicoes univariadas devem ser Normais ou LogNormais")
    else
        return true
    end
end

"""
    validate_univar_gaussian(models_per_stage::Vector{Vector{<:Distribution}})

Versao multivariada: garante que toda marginal, em todos os estagios e usinas, seja `Normal`
ou `LogNormal`.
"""
function validate_univar_gaussian(models_per_stage::Vector{Vector{<:Distribution}})
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

"""
    validate_copulas_kind(copulas_per_stage)

Garante que todas as copulas sejam `IndependentCopula` ou `GaussianCopula` — as unicas
compativeis com a construcao da incerteza gaussiana multivariada. Lanca erro caso contrario.
"""
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

"""
    validate_stch_for_ldr(stch::Naive)

Valida o processo estocastico para uso no LDR e retorna a tupla `(models, copulas)` com as
distribuicoes marginais e copulas por estagio ja verificadas por
[`validate_univar_gaussian`](@ref) e [`validate_copulas_kind`](@ref).
"""
function validate_stch_for_ldr(stch::Naive)
    models = get_models_per_stage(stch)
    copulas = get_copulas_per_stage(stch)
    #aux = validate_univar_gaussian(models) && validate_copulas_kind(copulas)
    return models, copulas
end

"""
    copula2Σ(c::IndependentCopula)

Converte uma copula independente na sua contribuicao para a matriz de covariancia: uma matriz
diagonal nula (sem correlacao cruzada entre as afluencias).
"""
function copula2Σ(c::IndependentCopula)
    diagm(repeat([0], length(c)))
end

"""
    copula2Σ(c::GaussianCopula)

Extrai da copula gaussiana a matriz de correlacao apenas com os termos fora da diagonal
(diagonal zerada), a ser somada aos desvios-padrao marginais na montagem da covariancia.

!!! note
    Modifica `c.Σ` in-place ao zerar a diagonal.
"""
function copula2Σ(c::GaussianCopula)
    Σ = c.Σ
    Σ[diagind(Σ)] .= 0.0
    return Σ
end

"""
    build_stage_ldruncertainty(stage_model::Distribution, stage_copula::Copula; kwargs...)

Constroi a distribuicao de incerteza univariada de um estagio, aplicando truncamento a
marginal. O truncamento e necessario porque o `Uncertainty` do LinearDecisionRules exige
suporte finito.

# Argumentos por chave
- `truncate_mode`: par `[modo_inf, modo_sup]`, cada um `:value` (limite absoluto) ou
  `:quantile` (quantil da propria distribuicao).
- `truncate_at`: par de valores/quantis correspondentes aos limites inferior e superior.
"""
function build_stage_ldruncertainty(
    stage_model::Distribution,
    stage_copula::Copula
    )::UnivariateDistribution

    return stage_model
end

"""
    build_stage_ldruncertainty(stage_models::Vector{<:Distribution}, stage_copula::IndependentCopula; kwargs...)

Constroi a incerteza de um estagio quando as afluencias sao independentes: trunca cada marginal
e as combina em uma `ProductDistribution`.
"""
function build_stage_ldruncertainty(
    stage_models::Vector{T} where T <: Distribution,
    stage_copula::IndependentCopula
    )::ProductDistribution

    return Distributions.product_distribution(stage_models...)
end

"""
    build_stage_ldruncertainty(stage_models::Vector{<:Distribution}, stage_copula::Copula)

Constroi a incerteza de um estagio com dependencia gaussiana: monta a matriz de covariancia a
partir dos desvios marginais somados a `copula2Σ(stage_copula)` e retorna uma `MvNormal` (se as
marginais forem `Normal`) ou `MvLogNormal` (se `LogNormal`).
"""
function build_stage_ldruncertainty(
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

"""
    stch2ldruncertainty(stch::Naive)::Vector{Distribution}

Converte o processo estocastico `Naive` do lab no vetor de distribuicoes de incerteza (uma por
estagio) esperado pelo LinearDecisionRules, validando-o e aplicando
[`build_stage_ldruncertainty`](@ref) a cada par (marginais, copula).
"""
function stch2ldruncertainty(stch::Naive)::Vector{Distribution}
    models, copulas = validate_stch_for_ldr(stch)

    mvdists = []
    for (m,c) in zip(models, copulas)
        push!(mvdists, build_stage_ldruncertainty(m, c))
    end

    return mvdists
end

"""
    determinisitic_stage1(dist::Distribution)::Vector{Float64}

Retorna o valor deterministico do primeiro estagio como a media da distribuicao, sempre na
forma de vetor (escalar embrulhado para distribuicoes univariadas). O primeiro estagio nao tem
incerteza no LDR, por isso a afluencia e fixada nesse valor esperado.
"""
function determinisitic_stage1(dist::Distribution)::Vector{Float64}
    isa(dist, UnivariateDistribution) ? [mean(dist)] : mean(dist)
end

"""
    add_inflow!(m::LDRModel, stch::Naive; n_breakpoints = 0)

Declara as variaveis de afluencia do modelo. O primeiro estagio recebe uma afluencia
deterministica, fixada na media (`determinisitic_stage1`); os demais recebem variaveis de
incerteza via [`add_uncertain_inflow!`](@ref).

Cada estagio possui sua propria variavel `INFLOW_s` (escalar quando ha uma unica usina, vetorial
caso contrario), pois nao ha como declarar uma unica variavel matricial associando cada coluna a
um conjunto de incerteza distinto.
"""
function add_inflow!(m::LDRModel, stch::Naive; n_breakpoints = 0)
    num_hydros, num_stages = size(stch)

    ldr_dists = stch2ldruncertainty(stch)
    stage_1 = determinisitic_stage1(ldr_dists[1])

    # ate entao nao tem como declarar uma unica INFLOW matricial colocando cala coluna num Uncertainty
    # set diferente
    # Em funcao disso, estagio precisa ter a propria variavel de inflow
    if num_hydros == 1
        m[:INFLOW_1] = JuMP.@variable(m, base_name = "INFLOW_1")
        JuMP.fix(m[:INFLOW_1], stage_1[1])
    else
        m[:INFLOW_1] = JuMP.@variable(m, [1:num_hydros], base_name = "INFLOW_1")
        for n in 1:num_hydros
            JuMP.fix(m[:INFLOW_1][n], stage_1[n])
        end
    end

    for s in 2:num_stages
        add_uncertain_inflow!(m, s, ldr_dists[s], n_breakpoints = n_breakpoints)
    end
end

"""
    add_uncertain_inflow!(m::LDRModel, stage::Int, dist::UnivariateDistribution; n_breakpoints = 0)

Declara a afluencia incerta escalar do estagio `stage` como variavel `INFLOW_stage` com o
conjunto `Uncertainty(distribution = dist)`. Quando `n_breakpoints > 0`, atribui pontos de
quebra (`BreakPoints`) para que a regra de decisao associada seja aproximada por partes.
"""
function add_uncertain_inflow!(m::LDRModel, stage::Int, dist::UnivariateDistribution;
    n_breakpoints = 0)

    varname = Symbol("INFLOW_" * string(stage))
    m[varname] = JuMP.@variable(
        m,
        base_name = String(varname),
        set = Uncertainty(; distribution = dist)
    )

    if n_breakpoints > 0
        JuMP.set_attribute(m[varname], BreakPoints(), n_breakpoints)
    end

    return nothing
end

"""
    add_uncertain_inflow!(m::LDRModel, stage::Int, dist::ProductDistribution; n_breakpoints = 0)

Declara a afluencia incerta vetorial (afluencias independentes) do estagio `stage` como um vetor
de variaveis `INFLOW_stage` sob o conjunto `Uncertainty`, opcionalmente com pontos de quebra por
componente.
"""
function add_uncertain_inflow!(m::LDRModel, stage::Int, dist::ProductDistribution;
    n_breakpoints = 0)

    size = length(dist)
    varname = Symbol("INFLOW_" * string(stage))
    m[varname] = JuMP.@variable(
        m,
        [1:size],
        base_name = String(varname),
        set = Uncertainty(; distribution = dist)
    )

    if n_breakpoints > 0
        for n in 1:size
            JuMP.set_attribute(m[varname][n], BreakPoints(), n_breakpoints)
        end
    end

    return nothing
end

"""
    add_uncertain_inflow!(m::LDRModel, stage::Int, dist::MultivariateNormal; n_breakpoints = 0, alpha = .95)

Declara a afluencia incerta vetorial com dependencia gaussiana do estagio `stage`. A `MvNormal`
e envolvida em uma `ConfidenceMvNormal` de nivel `alpha` — que a restringe ao elipsoide de
confianca, garantindo suporte finito exigido pelo `Uncertainty` — antes de declarar o vetor de
variaveis, opcionalmente com pontos de quebra por componente.
"""
function add_uncertain_inflow!(m::LDRModel, stage::Int, dist::MultivariateNormal;
    n_breakpoints = 0, alpha = .95)
    size = length(dist)
    dist = LinearDecisionRules.ConfidenceMvNormal(dist.μ, dist.Σ, alpha)
    varname = Symbol("INFLOW_" * string(stage))
    m[varname] = JuMP.@variable(
        m,
        [1:size],
        base_name = String(varname),
        set = Uncertainty(; distribution = dist)
    )

    if n_breakpoints > 0
        for n in 1:size
            JuMP.set_attribute(m[varname][n], BreakPoints(), n_breakpoints)
        end
    end

    return nothing
end

# END CUTS -----------------------------------------------------------------------------------------

function add_endcuts!(m::LDRModel, cuts::EndCuts)
    n_cuts = cuts.n_cuts

    m[:theta] = JuMP.@variable(m, base_name = "theta")
    JuMP.set_lower_bound(m[:theta], 0.0)

    JuMP.@constraint(
        m, [k in 1:n_cuts],
        m[:theta] >= cuts.intercept[k] + cuts.coefficients[k,:]' * m[STORED_VOLUME][:,end]
    )
end

# OBJECTIVE FUNCTION -------------------------------------------------------------------------------

"""
    add_objective!(m::LDRModel, s::SystemData, stch::Naive)

Define a funcao objetivo de minimizacao do custo total ao longo de todos os estagios, somando o
custo de geracao termica, a penalidade de deficit, as penalidades de intercambio direto e
reverso e a penalidade de vertimento.
"""
function add_objective!(m::LDRModel, s::SystemData, stch::Naive)
    hydros = get_hydros_entities(s)
    buses = get_buses_entities(s)
    lines = get_lines_entities(s)
    thermals = get_thermals_entities(s)
    nb = length(buses)
    nl = length(lines)
    nh = length(hydros)
    nt = length(thermals)
    ns = size(stch, 2)

    @objective(
        m,
        Min,
            m[:theta] +
            sum(thermals[n].cost * m[THERMAL_GENERATION][n,s] for n in 1:nt, s in 1:ns) +
            sum(buses[n].deficit_cost * m[DEFICIT][n,s] for n in 1:nb, s in 1:ns) +
            sum(lines[n].exchange_penalty * m[DIRECT_EXCHANGE][n,s] for n in 1:nl, s in 1:ns) +
            sum(lines[n].exchange_penalty * m[REVERSE_EXCHANGE][n,s] for n in 1:nl, s in 1:ns) +
            sum(hydros[n].spillage_penalty * m[SPILLAGE][n,s] for n in 1:nh, s in 1:ns)
    )
end
