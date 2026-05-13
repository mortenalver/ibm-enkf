using Random

mutable struct Individual
    x::Float64
    y::Float64
    v_x::Float64
    v_y::Float64
    n::Float64
end


function createIndividual(x, y, n)
    ind = Individual(x, y, 0.0, 0.0, 1)
    return ind
end

function step(t, dt, ind, indsArray, perturb, idx, xall, yall, X_fld, xrng, yrng, cog, ms)


end