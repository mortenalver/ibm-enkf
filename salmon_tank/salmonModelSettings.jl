
# Struct holding model settings:
mutable struct ModelSettings
    nInd::Int
    nPerInd::Float64
    radius::Float64
    xMax::Float64
    yMax::Float64
    dxy::Float64
    tau::Float64
    BL::Float64
    charVel::Float64
    d_wpref::Float64
    sigma::Float64
    d_l::Float64
    d_h::Float64
    ModelSettings() = new(
        80,
        1.0,
        8.0,
        8.0,
        4.0,
        0.5,
        0.9, # tau
        0.2, # BL, m
        0.4, # characteristic speed, BL/s
        0.25, # d_wpref (m)
        0.25,   # sigma
        0.66, # d_l (BL)
        3.0,  # d_h (BL)
        )
end