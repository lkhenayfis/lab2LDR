module lab2ldr

using LinearDecisionRules

import SDDPlab.Lab: Model, Engine
import SDDPlab.Engines: build

struct LDREngine <: Engine end

struct LDRModel <: Model
    model::LinearDecisionRules.LDRModel
end

end
