

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
    foodLimit::Float64
    fuzzySinkhornMoves::Bool
    anamorphicTransform::Bool
    AssimSettings() = new(false,2,10,false,true,false,false,1000,false,1,1.0,false,4.0,4.0,false,false)
end
