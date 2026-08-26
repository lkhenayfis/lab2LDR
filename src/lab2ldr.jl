module lab2ldr

using Distributions, Copulas, LinearAlgebra
using LinearDecisionRules
using JuMP

import LinearDecisionRules: LDRModel, Stage, Uncertainty, BreakPoints
import Copulas: Copula
import Distributions: ProductDistribution

# As of the abstract-engine branch, SDDPlab no longer re-exports these names from
# its top-level module; they live in (and are exported by) the submodules below.
# We import exactly the names this bridge uses.
import SDDPlab: read_study, Study
import SDDPlab.Utils: read_jsonc
import SDDPlab.System:
    Buses,
    Lines,
    Thermals,
    Hydros,
    SystemData,
    get_system,
    get_hydros,
    get_buses,
    get_lines,
    get_thermals,
    get_hydros_entities,
    get_buses_entities,
    get_lines_entities,
    get_thermals_entities,
    get_ids,
    downstream
import SDDPlab.Scenarios: get_scenarios, get_graph, get_number_of_stages, get_load
import SDDPlab.StochasticProcess: Naive
import SDDPlab.Lab:
    InputModule,
    LOAD,
    DEFICIT,
    DIRECT_EXCHANGE,
    REVERSE_EXCHANGE,
    NET_EXCHANGE,
    THERMAL_GENERATION,
    THERMAL_GENERATION_COST,
    STORED_VOLUME,
    TURBINED_FLOW,
    SPILLAGE,
    OUTFLOW,
    HYDRO_GENERATION,
    HYDRO_BALANCE,
    LOAD_BALANCE

include("input.jl")
include("engine-ldr.jl")

export build_model

end
