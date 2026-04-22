# This is the main file. To run the simulation in a new Julia REPL, first 
# evaluate this file, then call the main function from the command line.


# The following statements must be run once to install each library:
#import Pkg; Pkg.add("Plots")
#import Pkg; Pkg.add("PyPlot")
#import Pkg; Pkg.add("NCDatasets")
#import Pkg; Pkg.add("Statistics")
#import Pkg; Pkg.add("OptimalTransport")
#import Pkg; Pkg.add("DataAssim")

using Plots
using NCDatasets
using Statistics
using LinearAlgebra
using DelimitedFiles
using OptimalTransport
using DataAssim
using Random 

include("settings.jl")
include("ibmModel.jl")
include("computeDensityField.jl")
include("enKF.jl")
include("util.jl")
include("measurementModel.jl")
include("applyCorrectionsIBM.jl")
include("applyCorrectionsResample.jl")
include("applyCorrectionsSinkhorn.jl")
include("applyCorrectionsSinkhornParallel.jl")
include("ibmAssimilation.jl")

function readCurrentField(filename)
    ds = NCDataset(filename,"r")
    v_u = ds["u_velocity"]
    uval = v_u[:,:,1,1]
    v_v = ds["v_velocity"]
    vval = v_v[:,:,1,1]
    return uval, vval
end

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

storageDir = "C:/temp/"
if ~isdir(storageDir)
    storageDir = "./"
end

function main(setDryrun, setResample, recordingTwin)
    

    # Basic settings:
    simnamePrefix = "twintest2_" 
    useRecordedTwin = true
    recordedTwinPrefix = storageDir*"d_r2twin_resample_"
    dt = 0.2 # Time step
    t_end = 150.2 #99.8#119.8 # Simulation end time
    storageInterval = 1

    # Assimilation settings:
    as = AssimSettings()
    as.N = 300# 100 # Number of ensemble members.
    

    #recordingTwin = true # True to record new twin:
    if recordingTwin
        # Set random seed to get deterministic outcome. Note: Threads.nthreads must be set to 1!
        Random.seed!(246235)

        t_end = 120.0
        useRecordedTwin = false
        storageInterval = 1
        simnamePrefix = "r1twin"
        as.N = 1
    end
    plotTimeStep = 40
    initFoodLevel = 1.0
    
    # Simulation parameters:
    ms = ModelSettings()
    ms.xMax = 20
    ms.yMax = 15
    ms.k_X = 2.5 # Strength of motion towards food gradient.
    ms.indsInteraction = true # If true, individuals will be repulsed from each other at close distances (not very optimized, so makes model slower)
    ms.indsAlignThresh = 0.3 # Distance threshold for individual alignment.
    ms.indsRepulseThresh = 0.08 # Distance threshold for individual repulsion.
    ms.indsAttractThresh = 0.75 # Distance threshold for individual attraction.
    ms.n_wuv = 4 # Number of perturbations for u/v
    ms.n_wX = 2 # Number of perturbations for X
    ms.sigma_uv = 1.0 # 3.0 # Intensity of u/v random perturbations
    ms.sigma_X = 0.8 # Intensity of X random perturbations
    ms.d_wuv = 5.0
    ms.d_wX = 6.0
    ms.indsAlignStrength = 0.01 # Strength of individual interaction
    ms.indsRepulseStrength = 0.4#0.08 # Strength of individual interaction
    ms.indsAttractStrength = 0.002 # Strength of individual interaction
    ms.pullTowardsCOG = false # If true, individuals will be pulled towards the center of gravity of the population
    ms.pullTowardsCOGStrength = 0.75#0.3 # Strength of the pull towards COG if activated
    ms.speedUpdateRate = 0.4 # 0.6 # Multiplier for speed update - lower means more intertia in speed updates
    ms.nInd = 5000 #6000 # Number of individuals
    ms.nPerInd = 1.0 # individuals per super individual
    
    # More assimilation settings:
    as.dryRun = setDryrun # If true, the assimilation process will be run but changes will not be applied.
    as.resampleAll = setResample # True to use resampling strategy instead of sinkhorn/resize strategy
    as.assimInterval = 10 # Time steps between each assimilation procedure
    as.speedsInStateVec = true # If true, include mean speed components per grid cell in the state vector.
    as.foodInStateVec = true # If true, include food field in the state vector
    as.measureFood = true # If true, include measurements of food
    as.regularMeasurements = true # If true, place measurements regularly at a given measSpacing
    as.measSpacing = 3 # If regular measurements, sets the measurement spacing
    as.nmeas = 1200 #2*800 # Number of randomly distributed meaurements 
    as.measVar = 0.25^2.0 # Squared measurement standard deviation
    # Not relevant when using ESTFK: 
    as.perturbMeasurements = true # True to perturb measurement matrix bin analysis step
    as.localizationDist = 6.5 #3.0 # Localization distance
    as.fuzzySinkhornMoves = false # If true, all movements are made with random +- half a cell distance
    
    # Modify sim name according to run mode:
    simname = simnamePrefix*"_"
    if as.dryRun
        simname = "d_"*simname
    end
    if as.resampleAll
        simname = simname*"resample_"
    end

    # Define domain area:
    xlim = [0, ms.xMax]
    ylim = [0, ms.yMax]
    dxy = 0.5 # Grid resolution
    
    
    # Number of time steps:
    ntimes = round(Int, t_end/dt)

    nstoretimes = round(Int, ntimes / storageInterval)
    storeCount = 0
    storeXYE_twin = fill(0.0, 7, ms.nInd, nstoretimes)
    storeXYE_e1 = fill(0.0, 7, ms.nInd, nstoretimes)
    eFillval = 0.0

    println("dt=", dt, ", t_end=",t_end,", steps=",ntimes)

    # Ensemble setup:
    Ndim = as.N + 1

    # If we are to use a recorded twin, load it now:
    rTwinX = []
    rTwinY = []
    rTwinU = []
    rTwinV = []
    rTwinE = []
    rTwinN = []

    if useRecordedTwin
        rTwinX = readdlm(recordedTwinPrefix*"twinX.csv", ',')
        println(size(rTwinX))
        rTwinY = readdlm(recordedTwinPrefix*"twinY.csv", ',')
        println(size(rTwinY))
        rTwinVX = readdlm(recordedTwinPrefix*"twinVX.csv", ',')
        println(size(rTwinVX))
        rTwinVY = readdlm(recordedTwinPrefix*"twinVY.csv", ',')
        println(size(rTwinVY))
        rTwinE = readdlm(recordedTwinPrefix*"twinE.csv", ',')
        rTwinN = readdlm(recordedTwinPrefix*"twinN.csv", ',')
        rTwinXfld = readdlm(recordedTwinPrefix*"twinXfld.csv", ',')
    end


    # Initialize ensemble:
    ensemble = Array{Array{Individual,1},1}(undef, Ndim)
    for ensI = 1:Ndim
        ensemble[ensI] = Array{Individual, 1}(undef, ms.nInd)
    
        if ensI < Ndim || !useRecordedTwin         
            for i = 1:ms.nInd
                
                indX = 3.1 + 2*randn(Float64)
                indY = 9.1 + 2*randn(Float64)
                ensemble[ensI][i] = createIndividual(indX, indY, ms.nPerInd)
            end
        else
            # This is the twin and we should initialize from pre-recorded data:
            for i = 1:ms.nInd
                ensemble[ensI][i] = Individual(rTwinX[i,1], rTwinY[i,1], rTwinVX[i,1], rTwinVY[i,1], 
                    0.0, rTwinN[i,1], 0.0, 0.0, -1)
            end
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
    storeEnKF_field = fill(0.0, length(densityField), nstoretimes)
    storeStd_field = fill(0.0, length(densityField), nstoretimes)
    storeX_e = fill(0.0, length(densityField), nstoretimes)

    
    # Initialize food field on same dimensions as the density field:
    X_fld = fill(initFoodLevel, size(densityField,1), size(densityField,2),Ndim)
    
    # Let the twin's initial food field have a maximum:
    ##X_twin_pert = getRandomField([size(X_fld,1) size(X_fld,2)], 0.4, 1, size(X_fld,1)/2)
    #X_twin_pert = getNormalField([size(X_fld,1) size(X_fld,2)], 0.4, [15, 8], size(X_fld,1)/2)
    #X_fld[:,:,Ndim] += X_twin_pert

    # Initialize variable to hold the updated density field from the last
    # EnKF run (before application to the IBM)
    enkfField = zeros(Float64,size(densityField,1), size(densityField,2))

    # Main loop:
    for tstep = 1:ntimes

        t = (tstep-1)*dt
        println(0.1*round(Int, 10.0*t))

        # Time step of IBM:
        totInteractions = 0
            
        Threads.@threads for ensI = 1:Ndim
            if ensI < Ndim || !useRecordedTwin
                # Roll random numbers used to perturb individuals' speeds:
                perturb = zeros(Float64, ms.n_wuv, 4)
                for ptI = 1:size(perturb,1)
                    perturb[ptI,:] = [ms.sigma_uv*randn(Float64), ms.sigma_uv*randn(Float64),
                        xlim[1]+rand(Float64)*(xlim[2]-xlim[1]), ylim[1]+rand(Float64)*(ylim[2]-ylim[1])]
                end
                
                indsArray = ensemble[ensI]
                X_fld_upd, avgInteractions = stepAll(t, dt, indsArray, perturb, ms, X_fld[:,:,ensI], xrng, yrng)
                X_fld[:,:,ensI] = X_fld_upd
                totInteractions += avgInteractions
            else
                # This is the twin and we are using pre-recorded data:
                indsArray = ensemble[ensI] 
                for i = 1:length(indsArray)
                    ind = indsArray[i]
                    ind.x = rTwinX[i,tstep]
                    ind.y = rTwinY[i,tstep]
                    ind.v_x = rTwinVX[i,tstep]
                    ind.v_y = rTwinVY[i,tstep]
                    ind.n = rTwinN[i,tstep]
                    ind.E = rTwinE[i,tstep]
                end
                #X_fld_upd = stepAll(t, dt, indsArray, perturb, ms, X_fld[:,:,ensI], xrng, yrng)
                X_fld[:,:,ensI] = reshape(rTwinXfld[:,tstep], size(X_fld,1), size(X_fld,2))
            end

        end
        println("Average interactions per ind: "*string(totInteractions/as.N))

        if (mod(tstep, as.assimInterval) == 0) && !recordingTwin
            println("Assim at tstep=", tstep)
            
            doPlot = tstep==plotTimeStep
            updatedEnsemble, X_fld_upd, enkfField = ibmAssimilation(as, deepcopy(ensemble), copy(X_fld), xlim, ylim, dxy, doPlot)
            
            if !as.dryRun
                ensemble = updatedEnsemble
                X_fld = X_fld_upd
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
                # Find grid cell to store local food concentration:
                ix, iy = getGridCell(ind, xlim, ylim, dxy)
                if ix>=1 && iy>=1 && ix<=size(X_fld,1) && iy<=size(X_fld,2)
                    storeXYE_twin[5,indI,storeCount] = X_fld[ix, iy,Ndim]
                else
                    storeXYE_twin[5,indI,storeCount] = NaN
                end
                storeXYE_twin[6,indI,storeCount] = ind.v_x
                storeXYE_twin[7,indI,storeCount] = ind.v_y
            end
            indsArray = ensemble[1]
            for indI = 1:ms.nInd
                ind = indsArray[indI]
                storeXYE_e1[1,indI,storeCount] = ind.x
                storeXYE_e1[2,indI,storeCount] = ind.y
                storeXYE_e1[3,indI,storeCount] = ind.E
                storeXYE_e1[4,indI,storeCount] = ind.n
                # Find grid cell to store local food concentration:
                ix, iy = getGridCell(ind, xlim, ylim, dxy)
                if ix>=1 && iy>=1 && ix<=size(X_fld,1) && iy<=size(X_fld,2)
                    storeXYE_e1[5,indI,storeCount] = X_fld[ix, iy,1]
                else
                    storeXYE_e1[5,indI,storeCount] = NaN
                end
                storeXYE_e1[6,indI,storeCount] = ind.v_x
                storeXYE_e1[7,indI,storeCount] = ind.v_y
            end

            # Store density field for twin:
            densityField, xrng, yrng = computeDensityField(ensemble[Ndim], xlim, ylim, dxy)
            storeDens_twin[:,storeCount] = reshape(densityField, length(densityField), 1)
            energyField, xrng, yrng = computeAverageEnergyField(ensemble[Ndim], xlim, ylim, dxy, eFillval)
            storeEnergy_twin[:,storeCount] = reshape(energyField, length(densityField), 1)
            tmpU, tmpV = computeAverageSpeedField(ensemble[Ndim],xlim, ylim, dxy, 0.0)
            storeU_twin[:,storeCount] = reshape(tmpU, length(tmpU), 1)
            storeV_twin[:,storeCount] = reshape(tmpV, length(tmpV), 1)

            # Store latest EnKF field:
            storeEnKF_field[:,storeCount] = reshape(enkfField, length(enkfField),1)


            # Store X field for twin:
            storeX_twin[:,storeCount] = reshape(X_fld[:,:,Ndim], length(densityField), 1)
            # Store mean X field for ensemble:
            storeX_e[:,storeCount] = reshape(mean(X_fld[:,:,1:as.N],dims=3), length(densityField), 1)

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
            storeStd_field[:,storeCount] = std(densEnsemble, dims=2)
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

        # Flush stdout in case we are running on a HPC machine:
        flush(stdout)
    end
    
    # Storage directory:
    prefix = storageDir*simname
    
    # Store a single file giving the field dimensions and the assimilation interval divided by storage interval:
    writedlm(prefix*"fieldDims.csv", [size(densityField,1) size(densityField,2) as.assimInterval/storageInterval dt*storageInterval], ',')

    # Store twin states to files:
    writedlm(prefix*"twinX.csv", storeXYE_twin[1,:,:], ',')
    writedlm(prefix*"twinY.csv", storeXYE_twin[2,:,:], ',')
    writedlm(prefix*"twinVX.csv", storeXYE_twin[6,:,:], ',')
    writedlm(prefix*"twinVY.csv", storeXYE_twin[7,:,:], ',')
    writedlm(prefix*"twinE.csv", storeXYE_twin[3,:,:], ',')
    writedlm(prefix*"twinN.csv", storeXYE_twin[4,:,:], ',')
    writedlm(prefix*"twinFood.csv", storeXYE_twin[5,:,:], ',')
    writedlm(prefix*"twinDens.csv", storeDens_twin, ',')
    writedlm(prefix*"twinEnergy.csv", storeEnergy_twin, ',')
    writedlm(prefix*"twinU.csv", storeU_twin, ',')
    writedlm(prefix*"twinV.csv", storeV_twin, ',')
    writedlm(prefix*"twinXfld.csv", storeX_twin, ',')

    writedlm(prefix*"enkfField.csv", storeEnKF_field, ',')

    writedlm(prefix*"e1X.csv", storeXYE_e1[1,:,:], ',')
    writedlm(prefix*"e1Y.csv", storeXYE_e1[2,:,:], ',')
    writedlm(prefix*"e1E.csv", storeXYE_e1[3,:,:], ',')
    writedlm(prefix*"e1N.csv", storeXYE_e1[4,:,:], ',')
    writedlm(prefix*"e1Food.csv", storeXYE_e1[5,:,:], ',')
    writedlm(prefix*"e1VX.csv", storeXYE_e1[6,:,:], ',')
    writedlm(prefix*"e1VY.csv", storeXYE_e1[7,:,:], ',')

    writedlm(prefix*"eDens.csv", storeDens_e, ',')
    writedlm(prefix*"eDensStd.csv", storeStd_field, ',')
    writedlm(prefix*"eEnergy.csv", storeEnergy_e, ',')
    writedlm(prefix*"eXfld.csv", storeX_e, ',')
    
end



