

using Plots
using Statistics
using LinearAlgebra
using DelimitedFiles
using OptimalTransport
using DataAssim
using Random 

# Model code to include:

# Salmon model:
include("salmon_tank/salmonModel.jl")
include("salmon_tank/salmonModelSettings.jl")
include("salmon_tank/salmonMeasurementModel.jl")
simnamePrefix = "salm" 

# Test model:
#include("testModel.jl")
#include("testModelSettings.jl")
#include("testMeasurementModel.jl")
#simnamePrefix = "test" 


include("settings.jl")
include("computeDensityField.jl")
include("enKF.jl")
include("util.jl")

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

storageDir = "C:/temp/"
if ~isdir(storageDir)
    storageDir = "./"
end

totalTimeUsed = 0

function main(setDryrun, setResample, recordingTwin)

    useRecordedTwin = false
    recordedTwinPrefix = storageDir*"d_twin123_resample_"
    if recordingTwin
        useRecordedTwin = false
    end

    # Basic settings:
    dt = 0.2 # Time step
    t_end = 70.2#99.8#119.8 # Simulation end time
    storageInterval = 1

    plotTimeStep = 20


    # Simulation parameters:
    println(getModelName())
    ms = ModelSettings()
    setModelSettings(ms) # Let the model set its preferred settings.

    # Assimilation settings:
    as = AssimSettings()
    as.N = 3# 100 # Number of ensemble members.
    as.resampleAll = setResample
    as.dryRun = setDryrun
    as.foodInStateVec = true

    # Define domain area:
    xlim = [0, ms.xMax]
    ylim = [0, ms.yMax]
    dxy = ms.dxy # Grid resolution

    # Number of time steps:
    ntimes = round(Int, t_end/dt)

    nstoretimes = round(Int, ntimes / storageInterval)
    storeCount = 0
    println("dt=", dt, ", t_end=",t_end,", steps=",ntimes)
    println("xmax=", ms.xMax, ", ymax=", ms.yMax,", dxy=",dxy)
    # If we are to use a recorded twin, load it now:
    rTwinX = []
    rTwinY = []
    rTwinU = []
    rTwinV = []
    rTwinE = []
    rTwinN = []

    if useRecordedTwin
        rTwinX = readdlm(recordedTwinPrefix*"twinX.csv", ',')
        rTwinY = readdlm(recordedTwinPrefix*"twinY.csv", ',')
        rTwinVX = readdlm(recordedTwinPrefix*"twinVX.csv", ',')
        rTwinVY = readdlm(recordedTwinPrefix*"twinVY.csv", ',')
        rTwinE = readdlm(recordedTwinPrefix*"twinE.csv", ',')
        rTwinN = readdlm(recordedTwinPrefix*"twinN.csv", ',')
        rTwinXfld = readdlm(recordedTwinPrefix*"twinXfld.csv", ',')
    end

    # Initialize ensemble:
    Ndim = as.N + 1

    ensemble = Array{Array{Individual,1},1}(undef, Ndim)
    for ensI = 1:Ndim
        ensemble[ensI] = Array{Individual, 1}(undef, ms.nInd)
    
        if ensI < Ndim || !useRecordedTwin         
            for i = 1:ms.nInd
                
                startPos = getRandomStartPos(ms) # Query the model for a random starting position
                ensemble[ensI][i] = createIndividual(startPos[1], startPos[2], ms.nPerInd, ms)
            end
        else
            # This is the twin and we should initialize from pre-recorded data:
            for i = 1:ms.nInd
                ensemble[ensI][i] = createIndividual([rTwinX[i,1], rTwinY[i,1], rTwinVX[i,1], rTwinVY[i,1], 
                    0.0, rTwinN[i,1], rTwinE[i,1]])
            end
        end
    end

    # Compute an initial density field to get the dimensions:
    densityField, xrng, yrng = computeDensityField(ensemble[Ndim], xlim, ylim, dxy)

    # Initialize common field on same dimensions as the density field:
    X_fld = fill(getInitialFieldLevel(), size(densityField,1), size(densityField,2),Ndim)
    
    # Initialize fields for storage:
    storeDens_twin = fill(0.0, length(densityField), nstoretimes)
    storeU_twin = fill(0.0, length(densityField), nstoretimes)
    storeV_twin = fill(0.0, length(densityField), nstoretimes)
    storeEnergy_twin = fill(0.0, length(densityField), nstoretimes)
    storeX_twin = fill(0.0, length(densityField), nstoretimes)
    storeDens_e = fill(0.0, length(densityField), nstoretimes)
    storeStd_field = fill(0.0, length(densityField), nstoretimes)
    storeX_e = fill(0.0, length(densityField), nstoretimes)
    storeXYE_twin = fill(0.0, 7, ms.nInd, nstoretimes)
    storeXYE_e1 = fill(0.0, 7, ms.nInd, nstoretimes)
    storeEnergy_e = fill(0.0, length(densityField), nstoretimes)
    storeEnKF_field = fill(0.0, length(densityField), nstoretimes)
    eFillval = 0.0

    # Initialize variable to hold the updated density field from the last
    # EnKF run (before application to the IBM)
    enkfField = zeros(Float64,size(densityField,1), size(densityField,2))


    # Main loop:
    for tstep = 1:ntimes
        t = (tstep-1)*dt
        println(0.1*round(Int, 10.0*t))

        Threads.@threads for ensI = 1:Ndim
            if ensI < Ndim || !useRecordedTwin
                # Ordinary ensemble member, or non-prerecorded twin model:
                indsArray = ensemble[ensI]
                X_fld_upd, numval = stepAll(t, dt, indsArray, true, ms, X_fld[:,:,ensI], xrng, yrng)
                if !isnothing(X_fld_upd)
                    X_fld[:,:,ensI] = X_fld_upd
                end

                #println(size(X_fld))
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
                X_fld[:,:,ensI] = reshape(rTwinXfld[:,tstep], size(X_fld,1), size(X_fld,2))
            end

        end
    
        if (mod(tstep, as.assimInterval) == 0) && !recordingTwin
            println("Assim at tstep=", tstep)
            
            doPlot = tstep==plotTimeStep
            updatedEnsemble, X_fld_upd, enkfField = ibmAssimilation(as, deepcopy(ensemble), copy(X_fld), xlim, ylim, dxy, doPlot, t)
            
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

        end


        # Flush stdout in case we are running on a HPC machine:
        flush(stdout)
    end

    # Modify sim name according to run mode:
    simname = simnamePrefix*"_"
    if as.dryRun
        simname = "d_"*simname
    end
    if as.resampleAll && !as.dryRun
        simname = simname*"resample_"
    end

    prefix = storageDir*simname
    
    print(xlim)
    print(ylim)
    print(dxy)

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


