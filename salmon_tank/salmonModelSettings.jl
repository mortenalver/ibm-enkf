
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
    feedingStart::Float64
    feedingEnd::Float64
    feedingRate::Float64
    pelletWeight::Float64
    ModelSettings() = new(
        500, # nInd
        1.0, # nPerInd
        4.0, # radius
        8.0, # xMax
        8.0, # yMax
        0.25, # dxy
        0.9, # tau
        0.25, # BL, m
        0.4, # characteristic speed, BL/s
        0.5, # d_wpref (m)
        0.25,   # sigma
        0.66, # d_l (BL)
        3.0,  # d_h (BL)
        5.0, # Feeding start
        25.0, # Feeding end
        1.0,  # Feeding rate (g/s)
        0.05  # Pellet weight (g)
        )
end