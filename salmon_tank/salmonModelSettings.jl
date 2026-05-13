
# Struct holding model settings:
mutable struct ModelSettings
    nInd::Int
    nPerInd::Float64
    xMax::Float64
    yMax::Float64
    ModelSettings() = new(2000,1.0,10.0,15.0)
end