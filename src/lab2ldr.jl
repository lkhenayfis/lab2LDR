module lab2ldr

# loads with import so exported names don't come to the module
# currently LinearDecisionRules does not export any names but if this eventually changes it will
# create a conflict with our LDRModel
using LinearDecisionRules
import LinearDecisionRules: LDRModel, Stage
using JuMP
using SDDPlab

struct LDREngine <: Engine end

struct LDRModel <: Model
    model::LinearDecisionRules.LDRModel
end

end
