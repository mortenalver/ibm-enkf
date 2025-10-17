# This is the main file. To run the simulation in a new Julia REPL, first 
# evaluate this file, then call the main function from the command line.


# The following statements must be run once to install each library:
#import Pkg; Pkg.add("Plots")
#import Pkg; Pkg.add("PyPlot")
#import Pkg; Pkg.add("NCDatasets")
#import Pkg; Pkg.add("Statistics")
#import Pkg; Pkg.add("OptimalTransport")

using Plots
using NCDatasets
using Statistics
using LinearAlgebra
using DelimitedFiles
using OptimalTransport

include("ibmModel.jl")
include("computeDensityField.jl")
include("enKF.jl")
include("util.jl")
include("measurementModel.jl")
include("applyCorrectionsIBM.jl")
include("applyCorrectionsResample.jl")
include("applyCorrectionsSinkhorn.jl")
include("ibmAssimilation.jl")

function readCurrentField(filename)
    ds = NCDataset(filename,"r")
    v_u = ds["u_velocity"]
    uval = v_u[:,:,1,1]
    v_v = ds["v_velocity"]
    vval = v_v[:,:,1,1]
    return uval, vval
end

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

    AssimSettings() = new(false,2,10,false,true)
end

function main(setDryrun, setResample)

    # Basic settings:
    simnamePrefix = "r13"
    dt = 0.1 # Time step
    t_end = 105.6 # Simulation end time
    storageInterval = 2
    initFoodLevel = 1.0
    
    # Simulation parameters:
    ms = ModelSettings()
    ms.migration = false # If true, periodic migration through the four corners of the domain. If false, migration controlled by food field and random movement.
    ms.indsInteraction = false # If true, individuals will be repulsed from each other at close distances (not very optimized, so makes model slower)
    ms.indsInteractionThresh = 0.22 # Distance threshold for individual interaction.
    ms.indsInteractionStrength = 0.1 # Strength of individual interaction
    ms.speedUpdateRate = 0.6 # Multiplier for speed update - lower means more intertia in speed updates
    ms.nInd = 2000 #6000 # Number of individuals
    ms.nPerInd = 1.0 # individuals per super individual
    ms.minNormSpeed = 0.5 # In migration mode, determines minimum typical speed of individuals.
    ms.scopeNormSpeed = 1.8 # In migration mode, determines scope of the typical speed of individuals.

    # Assimilation settings:
    as = AssimSettings()
    as.dryRun = setDryrun # If true, the assimilation process will be run but changes will not be applied.
    as.N = 100# 100 # Number of ensemble members.
    as.resampleAll = setResample # True to use resampling strategy instead of sinkhorn/resize strategy
    as.assimInterval = 15 # Time steps between each assimilation procedure
    as.speedsInStateVec = false # If true, include mean speed components per grid cell in the state vector.

    # Modify sim name according to run mode:
    simname = simnamePrefix*"_"
    if as.dryRun
        simname = "d_"*simname
    end
    if as.resampleAll
        simname = simname*"resample_"
    end

    # Define domain area:
    xlim = [0, 20]
    ylim = [0, 15]
    dxy = 0.5 # Grid resolution
    
    
    # Number of time steps:
    ntimes = round(Int, t_end/dt)

    nstoretimes = round(Int, ntimes / storageInterval)
    storeCount = 0
    storeXYE_twin = fill(0.0, 4, ms.nInd, nstoretimes)
    storeXYE_e1 = fill(0.0, 4, ms.nInd, nstoretimes)
    eFillval = 0.0

    println("dt=", dt, ", t_end=",t_end,", steps=",ntimes)

    # Ensemble setup:
    Ndim = as.N + 1

    # Initialize ensemble:
    ensemble = Array{Array{Individual,1},1}(undef, Ndim)
    for ensI = 1:Ndim
        ensemble[ensI] = Array{Individual, 1}(undef, ms.nInd)
    
        for i = 1:ms.nInd
            normspeed = ms.minNormSpeed+ms.scopeNormSpeed*rand(Float64)
            # Modify speed for twin to create an offset (only affects simulation with migration activated):
            if ensI != Ndim
                normspeed = 1.25*normspeed # Non-twin is 25% faster
            end

            indX = 3.5 + 1.5*randn(Float64)
            indY = 9.5 + 1.5*randn(Float64)
            ensemble[ensI][i] = createIndividual(indX, indY, normspeed, ms.nPerInd)
        end
    end

    # Compute an initial density field to get the dimensions:
    densityField, xrng, yrng = computeDensityField(ensemble[Ndim], xlim, ylim, dxy)
    storeDens_twin = fill(0.0, length(densityField), nstoretimes)
    storeEnergy_twin = fill(0.0, length(densityField), nstoretimes)
    storeU_twin = fill(0.0, length(densityField), nstoretimes)
    storeV_twin = fill(0.0, length(densityField), nstoretimes)
    storeX_twin = fill(0.0, length(densityField), nstoretimes)
    storeDens_e = fill(0.0, length(densityField), nstoretimes)
    storeEnergy_e = fill(0.0, length(densityField), nstoretimes)

    # Initialize food field on same dimensions as the density field:
    X_fld = fill(initFoodLevel, size(densityField,1), size(densityField,2),Ndim)


    # Main loop:
    for tstep = 1:ntimes

        t = (tstep-1)*dt
        println(0.1*round(Int, 10.0*t))

        # Time step of IBM:
        for ensI = 1:Ndim

            # Roll random numbers used to perturb individuals' speeds:
            perturb = zeros(Float64, 20, 4)
            for ptI = 1:size(perturb,1)
                perturb[ptI,:] = [1.0*randn(Float64), 1.0*randn(Float64),
                    xlim[1]+rand(Float64)*(xlim[2]-xlim[1]), ylim[1]+rand(Float64)*(ylim[2]-ylim[1])]
            end
            
            indsArray = ensemble[ensI]
            X_fld_upd = stepAll(t, dt, indsArray, perturb, ms, X_fld[:,:,ensI], xrng, yrng)
            X_fld[:,:,ensI] = X_fld_upd
        end

        if mod(tstep, as.assimInterval) == 0
            println("Assim at tstep=", tstep)
            
            doPlot = tstep==15
            updatedEnsemble = ibmAssimilation(as, deepcopy(ensemble), xlim, ylim, dxy, doPlot)
            
            if !as.dryRun
                ensemble = updatedEnsemble
            end
        end

        # Store the twin model at regular time steps to a temporary array:
        if mod(tstep, storageInterval) == 0
            storeCount = storeCount+1
            # Store individuals:
            indsArray = ensemble[Ndim]
            for indI = 1:ms.nInd
                ind = indsArray[indI]
                storeXYE_twin[1,indI,storeCount] = ind.x
                storeXYE_twin[2,indI,storeCount] = ind.y
                storeXYE_twin[3,indI,storeCount] = ind.E
                storeXYE_twin[4,indI,storeCount] = ind.n
            end
            indsArray = ensemble[1]
            for indI = 1:ms.nInd
                ind = indsArray[indI]
                storeXYE_e1[1,indI,storeCount] = ind.x
                storeXYE_e1[2,indI,storeCount] = ind.y
                storeXYE_e1[3,indI,storeCount] = ind.E
                storeXYE_e1[4,indI,storeCount] = ind.n
            end

            # Store density field for twin:
            densityField, xrng, yrng = computeDensityField(ensemble[Ndim], xlim, ylim, dxy)
            storeDens_twin[:,storeCount] = reshape(densityField, length(densityField), 1)
            energyField, xrng, yrng = computeAverageEnergyField(ensemble[Ndim], xlim, ylim, dxy, eFillval)
            storeEnergy_twin[:,storeCount] = reshape(energyField, length(densityField), 1)
            tmpU, tmpV = computeAverageSpeedField(ensemble[Ndim],xlim, ylim, dxy, 0.0)
            storeU_twin[:,storeCount] = reshape(tmpU, length(tmpU), 1)
            storeV_twin[:,storeCount] = reshape(tmpV, length(tmpV), 1)

            # Store X field for twin:
            storeX_twin[:,storeCount] = reshape(X_fld[:,:,Ndim], length(densityField), 1)


            # Store mean ensemble density field:
            densEnsemble = zeros(Float64,length(xrng)*length(yrng), as.N)
            energyEnsemble = zeros(Float64,length(xrng)*length(yrng), as.N)
            for ensI = 1:as.N
                indsArray = ensemble[ensI]
                densityField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)
                densEnsemble[:,ensI] = reshape(densityField, length(densityField), 1)
                
                energyField, xrng, yrng = computeAverageEnergyField(indsArray, xlim, ylim, dxy, eFillval)
                energyEnsemble[:,ensI] = reshape(energyField, length(densityField), 1)
            end    
            storeDens_e[:,storeCount] = mean(densEnsemble, dims=2)
            storeEnergy_e[:,storeCount] = mean(energyEnsemble, dims=2)

            #densityField, xrng, yrng = computeDensityField(ensemble[1], xlim, ylim, dxy)
            #energyField, xrng, yrng = computeAverageEnergyField(ensemble[1], xlim, ylim, dxy, eFillval)
            #storeDens_e[:,storeCount] = reshape(densityField, length(densityField), 1)
            #storeEnergy_e[:,storeCount] = reshape(energyField, length(densityField), 1)

        end

        # if tstep==110 || tstep==120 || tstep==130
            
        #     xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
        #     yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
        #     dimensions = [length(xrng), length(yrng)]
        #     densEnsemble = zeros(Float64,length(xrng)*length(yrng), N)
        #     for ensI = 1:N
        #         indsArray = ensemble[ensI]
        #         densityField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)
        #         densEnsemble[:,ensI] = reshape(densityField, length(densityField), 1)
        #     end    
        #     updMeanField = reshape(mean(densEnsemble, dims=2), dimensions[1], dimensions[2])
        #     display(heatmap(updMeanField))
        # end
    end
    
    # Storage directory:
    prefix = "C:/temp/"*simname
    
    # Store a single file giving the field dimensions:
    writedlm(prefix*"fieldDims.csv", size(densityField), ',')

    # Store twin states to files:
    writedlm(prefix*"twinX.csv", storeXYE_twin[1,:,:], ',')
    writedlm(prefix*"twinY.csv", storeXYE_twin[2,:,:], ',')
    writedlm(prefix*"twinE.csv", storeXYE_twin[3,:,:], ',')
    writedlm(prefix*"twinN.csv", storeXYE_twin[4,:,:], ',')
    writedlm(prefix*"twinDens.csv", storeDens_twin, ',')
    writedlm(prefix*"twinEnergy.csv", storeEnergy_twin, ',')
    writedlm(prefix*"twinU.csv", storeU_twin, ',')
    writedlm(prefix*"twinV.csv", storeV_twin, ',')
    writedlm(prefix*"twinXfld.csv", storeX_twin, ',')

    writedlm(prefix*"e1X.csv", storeXYE_e1[1,:,:], ',')
    writedlm(prefix*"e1Y.csv", storeXYE_e1[2,:,:], ',')
    writedlm(prefix*"e1E.csv", storeXYE_e1[3,:,:], ',')
    writedlm(prefix*"e1N.csv", storeXYE_e1[4,:,:], ',')
    writedlm(prefix*"eDens.csv", storeDens_e, ',')
    writedlm(prefix*"eEnergy.csv", storeEnergy_e, ',')
end



