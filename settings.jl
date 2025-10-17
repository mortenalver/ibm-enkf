
# Struct holding model settings:
mutable struct ModelSettings
    nInd::Int
    nPerInd::Float64
    speedUpdateRate::Float64
    indsInteraction::Bool
    migration::Bool
    indsInteractionThresh::Float64
    indsInteractionStrength::Float64
    minNormSpeed::Float64
    scopeNormSpeed::Float64

    ModelSettings() = new(2000,1,0.6,true,false,1.0,2,3)
end

# Struct holding assimilation settings:
mutable struct AssimSettings
    dryRun::Bool
    N::Int
    assimInterval::Int
    resampleAll::Bool
    speedsInStateVec::Bool
    nmeas::Int 

    AssimSettings() = new(false,2,10,false,true,1000)
end