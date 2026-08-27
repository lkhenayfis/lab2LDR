# CUTS ---------------------------------------------------------------------------------------------

"""
    EndCuts

Conjunto de cortes lineares (planos de corte) que representam a funcao de custo futuro sobre o
volume armazenado no ultimo estagio — a condicao de contorno terminal do problema.

Cada corte `k` tem a forma `θ >= intercept[k] + coefficients[k, :] ⋅ v_final`, onde `v_final` e
o vetor de volumes finais das usinas. Como o formato de entrada atual do lab ainda nao modela
cortes/condicao de contorno nativamente, eles sao lidos de um bloco `cuts` adicionado ao
`main.jsonc` do estudo.

# Campos
- `intercept::Vector{Float64}`: interceptos dos cortes (um por corte).
- `coefficients::Matrix{Float64}`: coeficientes dos cortes, dimensao `n_cuts × n_hidros`.
- `n_cuts::Int`: numero de cortes.
"""
struct EndCuts
    intercept::Vector{Float64}
    coefficients::Matrix{Float64}
    n_cuts::Int
end

"""
    EndCuts(intercept::Vector{Float64}, coefficients::Matrix{Float64})

Constroi um [`EndCuts`](@ref) inferindo `n_cuts` a partir do comprimento de `intercept`.
"""
function EndCuts(intercept::Vector{Float64}, coefficients::Matrix{Float64})
    n = length(intercept)
    EndCuts(intercept, coefficients, n)
end

"""
    EndCuts(path::String)

Le o bloco `cuts` do `main.jsonc` do estudo em `path` e constroi o [`EndCuts`](@ref)
correspondente. Espera as chaves `intercept` (vetor) e `coefficients` (lista de linhas, empilhada
em matriz `n_cuts × n_hidros` via `stack`).
"""
function EndCuts(path::String)
    arq = path * "/main.jsonc"
    d = read_jsonc(arq, CompositeException())

    cuts = d["cuts"]
    cuts["intercept"] = real(cuts["intercept"])
    cuts["coefficients"] = real(stack(cuts["coefficients"]))

    EndCuts(cuts["intercept"], cuts["coefficients"])
end

# TRUNCATION ---------------------------------------------------------------------------------------

"""
    get_truncation_bounds(v::Vector, id::Int, stage::Int)

Extrai o par `[inferior, superior]` de truncamento da afluencia da usina `id` no estagio `stage`
a partir da lista `v` do bloco `truncation` (cada elemento tem `id` e uma lista de `bounds` por
`stage`).
"""
function get_truncation_bounds(v::Vector, id::Int, stage::Int)
    elem_config = filter(x -> x["id"] == id, v)[1]
    elem_stage_bounds = filter(x -> x["stage"] == stage, elem_config["bounds"])[1]
    elem_stage_bounds["value"]
end

"""
    do_truncation!(path::String, study::Study)

Aplica in-place o truncamento das distribuicoes marginais de afluencia do processo estocastico
de `study`, conforme o bloco `truncation` do `main.jsonc` em `path`.

Para cada usina e estagio estocastico (`s` de 2 ate o numero de estagios), substitui a marginal
`stch.models[i].distributions[s]` pela sua versao `truncated(dist, bounds...)`. O truncamento e
necessario porque o `Uncertainty` do LinearDecisionRules exige suporte finito; como o formato de
entrada do lab ainda nao carrega esses limites, eles vem do bloco `truncation`.

!!! note
    Modifica o `study` recebido (as distribuicoes do `stochastic_process`).
"""
function do_truncation!(path::String, study::Study)
    arq = path * "/main.jsonc"
    d = read_jsonc(arq, CompositeException())
    stch = study.inputs.files[1].inflow.stochastic_process

    trunc = d["truncation"]

    num_elems, num_stages = size(stch)

    for i in 1:num_elems, s in 2:num_stages
        bounds = get_truncation_bounds(trunc, i, s)
        dist = stch.models[i].distributions[s]
        stch.models[i].distributions[s] = truncated(dist, bounds...)
    end
end
