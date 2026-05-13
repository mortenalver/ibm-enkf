
# Struct holding model settings:
mutable struct ModelSettings
    nInd::Int
    nPerInd::Float64
    speedUpdateRate::Float64
    k_X::Float64
    X_a::Float64
    indsInteraction::Bool
    indsAlignThresh::Float64
    indsRepulseThresh::Float64
    indsAttractThresh::Float64
    indsRepulseStrength::Float64
    indsAlignStrength::Float64
    indsAttractStrength::Float64
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
    ModelSettings() = new(2000,1,0.6,1.0,0.01,true,1.0,1.0,1.0,1.0,1.0,
        false,1.0,5,5,1.0,1.0,2.0,2.0,10.0,10.0)
end