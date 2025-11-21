
# Struct holding model settings:
mutable struct ModelSettings
    nInd::Int
    nPerInd::Float64
    speedUpdateRate::Float64
    k_X::Float64
    indsInteraction::Bool
    indsAlignThresh::Float64
    indsRepulseThresh::Float64
    indsRepulseStrength::Float64
    indsAlignStrength::Float64
    pullTowardsCOG::Bool
    pullTowardsCOGStrength::Float64
    n_wuv::Int
    n_wX::Int
    sigma_uv::Float64
    sigma_X::Float64
    d_wuv::Float64
    d_wX::Float64
    xMax::Float64
    yMax::Float64
    ModelSettings() = new(2000,1,0.6,1.0,true,1.0,1.0,1.0,
        false,1.0,5,5,1.0,1.0,2.0,2.0,10.0,10.0)
end

# Struct holding assimilation settings:
mutable struct AssimSettings
    dryRun::Bool
    N::Int
    assimInterval::Int
    resampleAll::Bool
    speedsInStateVec::Bool
    foodInStateVec::Bool
    measureFood::Bool
    nmeas::Int
    regularMeasurements::Bool
    measSpacing::Int
    measVar::Float64
    perturbMeasurements::Bool
    localizationDist::Float64

    AssimSettings() = new(false,2,10,false,true,false,false,1000,false,1,1.0,false,4.0)
end