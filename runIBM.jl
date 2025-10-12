#import Pkg; Pkg.add("Plots")
#import Pkg; Pkg.add("PyPlot")
#import Pkg; Pkg.add("NCDatasets")
#import Pkg; Pkg.add("Statistics")
using Plots
#using PyPlot
using NCDatasets
using Statistics
using LinearAlgebra
using DelimitedFiles

include("ibmModel.jl")
#include("ibmModel2.jl")
include("computeDensityField.jl")
include("enKF.jl")
include("util.jl")
include("measurementModel.jl")
include("applyCorrectionsIBM.jl")
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
    indsInteraction::Bool
    indsInteractionThresh::Float64
    minNormSpeed::Float64
    scopeNormSpeed::Float64

    ModelSettings() = new(2000,1,false,1.0,2,3)
end

# Struct holding assimilation settings:
mutable struct AssimSettings
    dryRun::Bool
    N::Int
    assimInterval::Int

    AssimSettings() = new(false,2,10)
end

function main()

    # Basic settings:
    simname = "r1_"
    dt = 0.1 # Time step
    t_end = 20 # Simulation end time
    storageInterval = 2

    # Simulation parameters:
    ms = ModelSettings()
    ms.indsInteraction = true
    ms.nInd = 2000 #6000 # Number of individuals
    ms.nPerInd = 1.0 # individuals per super individual
    ms.minNormSpeed = 0.5
    ms.scopeNormSpeed = 1.8

    # Assimilation settings:
    as = AssimSettings()
    as.dryRun = true
    as.N = 20# 100
    as.assimInterval = 15

    # Define area for density fields:
    xlim = [0, 20]
    ylim = [0, 15]
    dxy = 0.5#(xlim[2]-xlim[1])/100
    
    

    # Number of time steps:
    ntimes = round(Int, t_end/dt)

    nstoretimes = round(Int, ntimes / storageInterval)
    storeCount = 0
    storeXYE_twin = fill(0.0, 4, ms.nInd, nstoretimes)
    storeXYE_e1 = fill(0.0, 4, ms.nInd, nstoretimes)
    eFillval = 0.0#NaN
    #allX = zeros(nInd, ntimes)
    #allY = zeros(nInd, ntimes)
    println("dt=", dt, ", t_end=",t_end,", steps=",ntimes)

    # Ensemble setup:
    Ndim = as.N + 1

    # Initialize ensemble:
    ensemble = Array{Array{Individual,1},1}(undef, Ndim)
    for ensI = 1:Ndim
        ensemble[ensI] = Array{Individual, 1}(undef, ms.nInd)
    
        for i = 1:ms.nInd
            normspeed = ms.minNormSpeed+ms.scopeNormSpeed*rand(Float64)
            # Modify speed for twin to create an offset:
            if ensI != Ndim
                normspeed = 1.25*normspeed # Non-twin is 25% faster
            end

            ensemble[ensI][i] = createIndividual(.5 + 5*(rand(Float64)-0.5), 
                7.5 + 15*(rand(Float64)-0.5), normspeed, ms.nPerInd)
        end
    end

    # Compute an initial density field to get the dimensions:
    densityField, xrng, yrng = computeDensityField(ensemble[Ndim], xlim, ylim, dxy)
    storeDens_twin = fill(0.0, length(densityField), nstoretimes)
    storeEnergy_twin = fill(0.0, length(densityField), nstoretimes)
    storeDens_e = fill(0.0, length(densityField), nstoretimes)
    storeEnergy_e = fill(0.0, length(densityField), nstoretimes)

    for tstep = 1:ntimes

        t = (tstep-1)*dt
        println(0.1*round(Int, 10.0*t))
        # Time step of IBM:
        for ensI = 1:Ndim
            
            
            perturb = zeros(Float64, 20, 4)
            #if ensI < Ndim
            for ptI = 1:size(perturb,1)
                perturb[ptI,:] = [2.5*(rand(Float64)-0.5), 2.5*(rand(Float64)-0.5),
                    xlim[1]+rand(Float64)*(xlim[2]-xlim[1]), ylim[1]+rand(Float64)*(ylim[2]-ylim[1])]
            end
            #end

            indsArray = ensemble[ensI]
            stepAll(t, dt, indsArray, perturb, ms)

        end

        if mod(tstep, as.assimInterval) == 0
            println("Assim at tstep=", tstep)
            
            doPlot = tstep==60
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

    writedlm(prefix*"e1X.csv", storeXYE_e1[1,:,:], ',')
    writedlm(prefix*"e1Y.csv", storeXYE_e1[2,:,:], ',')
    writedlm(prefix*"e1E.csv", storeXYE_e1[3,:,:], ',')
    writedlm(prefix*"e1N.csv", storeXYE_e1[4,:,:], ',')
    writedlm(prefix*"eDens.csv", storeDens_e, ',')
    writedlm(prefix*"eEnergy.csv", storeEnergy_e, ',')
end



