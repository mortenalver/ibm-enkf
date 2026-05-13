
# Model code to include:

# Salmon model:
#include("salmon_tank/salmonModel.jl")
#include("salmon_tank/salmonModelSettings.jl")

# Test model:
include("ibmModel.jl")
include("testModelSettings.jl")

include("computeDensityField.jl")
include("enKF.jl")
include("util.jl")
#include("measurementModel.jl")
include("applyCorrectionsIBM.jl")
include("applyCorrectionsResample.jl")
include("applyCorrectionsSinkhorn.jl")
include("applyCorrectionsSinkhornParallel.jl")
include("ibmAssimilation.jl")


function dryrun()
    main(true, true, false)
end

function normrun()
    main(false, false, false)
end

function resamplerun()
    main(false, true, false)
end

function recordtwin()
    main(true, true, true)
end

function allrun()
    dryrun()
    resamplerun()
    normrun()
end


function main(setDryrun, setResample, recordingTwin)

    # Basic settings:
    simnamePrefix = "salm" 
    dt = 0.2 # Time step
    t_end = 35.2#99.8#119.8 # Simulation end time
    storageInterval = 1

    # Simulation parameters:
    ms = ModelSettings()
    println(ms)

    ## Assimilation settings:
    #as = AssimSettings()
    #as.N = 300# 100 # Number of ensemble members.

    # Define domain area:
    xlim = [0, ms.xMax]
    ylim = [0, ms.yMax]
    dxy = 0.5 # Grid resolution

    ind = createIndividual(5, 7, 1)
    println(ind)
end


