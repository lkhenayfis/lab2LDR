# CUTS ---------------------------------------------------------------------------------------------

struct EndCuts
    intercept::Vector{Float64}
    coefficients::Matrix{Float64}
    n_cuts::Int
end

function EndCuts(intercept::Vector{Float64}, coefficients::Matrix{Float64})
    n = length(intercept)
    EndCuts(intercept, coefficients, n)
end

function EndCuts(path::String)
    arq = path * "/main.jsonc"
    d = read_jsonc(arq, CompositeException())

    cuts = d["cuts"]
    cuts["intercept"] = real(cuts["intercept"])
    cuts["coefficients"] = real(stack(cuts["coefficients"]))

    EndCuts(cuts["intercept"], cuts["coefficients"])
end

# TRUNCATION ---------------------------------------------------------------------------------------

function get_truncation_bounds(v::Vector, id::Int, stage::Int)
    elem_config = filter(x -> x["id"] == id, v)[1]
    elem_stage_bounds = filter(x -> x["stage"] == stage, elem_config["bounds"])[1]
    elem_stage_bounds["value"]
end

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
