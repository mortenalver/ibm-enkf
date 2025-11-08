
function ibmAssimilation(as, ensemble, X_fld, xlim, ylim, dxy, doPlot)
    # Handle the assimilation process from derivation of state values
    # via calling the Ensemble Kalman filter to updating the IBM.

    println("doplot="*string(doPlot))

    Ndim = as.N + 1 # There are N ensemble members plus a twin

    # Define ranges:
    xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
    yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
    dimensions = [length(xrng), length(yrng)]
    nPos = length(xrng)*length(yrng)
    statesPerPos = 1
    stateOffsetSpeed = 1
    stateOffsetFood = -1
    if as.speedsInStateVec
        stateOffsetSpeed = statesPerPos
        statesPerPos += 2
    end
    if as.foodInStateVec
        stateOffsetFood = statesPerPos
        statesPerPos += 1
    end
    println("statesPerPos: "*string(statesPerPos))
    
    densEnsemble = zeros(Float64,statesPerPos*nPos, as.N)
    densTwin = zeros(Float64, statesPerPos*nPos, 1)

    println("size densTwin: "*string(size(densTwin)))

    # Set up x and y values per state for localization:
    x_states = zeros(Float64, statesPerPos*nPos, 1)
    y_states = zeros(Float64, statesPerPos*nPos, 1)
    cind = CartesianIndices((length(xrng), length(yrng)))
    for i=1:nPos
        xi = cind[i][1]
        yi = cind[i][2]
        x_states[i] = dxy*xi
        y_states[i] = dxy*yi
        if as.speedsInStateVec
            x_states[i+nPos] = dxy*xi
            y_states[i+nPos] = dxy*yi
            x_states[i+2*nPos] = dxy*xi
            y_states[i+2*nPos] = dxy*yi                
        end
    end
    
    for ensI = 1:Ndim
        indsArray = ensemble[ensI]
        densityField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)

        if ensI <= as.N
            densEnsemble[1:nPos,ensI] = reshape(densityField, length(densityField), 1)
        else
            densTwin[1:nPos] = reshape(densityField, length(densityField), 1)
        end

        # Add speeds to state vector if we are supposed to:
        if as.speedsInStateVec
            tmpU, tmpV = computeAverageSpeedField(indsArray, xlim, ylim, dxy, 0.0)
        
            if doPlot && ensI==1
                writedlm("C:/temp/Ubefore.csv", tmpU, ',')
                writedlm("C:/temp/Vbefore.csv", tmpV, ',')
            end

            if ensI <= as.N
                densEnsemble[stateOffsetSpeed*nPos+1:(stateOffsetSpeed+1)*nPos,ensI] = reshape(tmpU, length(densityField), 1)
                densEnsemble[(stateOffsetSpeed+1)*nPos+1:(stateOffsetSpeed+2)*nPos,ensI] = reshape(tmpV, length(densityField), 1)
            else
                densTwin[stateOffsetSpeed*nPos+1:(stateOffsetSpeed+1)*nPos] = reshape(tmpU, length(densityField), 1)
                densTwin[(stateOffsetSpeed+1)*nPos+1:(stateOffsetSpeed+2)*nPos] = reshape(tmpV, length(densityField), 1)
            end
        end

        if as.foodInStateVec && doPlot && ensI==1
            # Store ens member 1 food before:
            writedlm("C:/temp/Xbefore.csv", X_fld[:,:,1], ',')
        end

        # Add food to state vector if we are supposed to:
        if as.foodInStateVec
            if ensI <= as.N
                densEnsemble[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos,ensI] = reshape(X_fld[:,:,ensI], length(densityField), 1)
            else
                densTwin[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos] = reshape(X_fld[:,:,ensI], length(densityField), 1)
            end
        end
    end    

    println("Member 1 sum abs deviation: ", sum(abs.(densEnsemble[:,1]-densTwin)) )

    # Ensemble Kalman filter:

    
    # Measurement model:
    M = zeros(Float64,1,1)
    xmeas = zeros(Float64,1,1)
    ymeas = zeros(Float64,1,1)
    if as.regularMeasurements
        M, xmeas, ymeas = getMRegular(dimensions, dxy, as, statesPerPos)
    else
        M = getM(dimensions, as, statesPerPos) 
    end
    xloc = getLocMatrix(dimensions, M, as.localizationDist)

    writedlm("C:/temp/Xloc.csv", xloc, ',')
    writedlm("C:/temp/M.csv", M, ',')

    y = M*densTwin # Measurement vector based on twin
    Rval = as.measVar # Assumed measurement uncertainty
    
    Rvec = Rval*ones(Float64, size(M,1))
    R = Diagonal(Rvec)
    


    #X_upd, X_upd_mean = DataAssim.ESTKF(densEnsemble, M*densEnsemble, y, R, M)

    #selectObs(i) = compact_locfun(sqrt.((x_states[i] .- xmeas).^2 .+ (y_states[i] .- ymeas).^2))#./as.localizationDist)
    #selectObs(i) = Vector(ones(Float64, length(xmeas))) #exp.(-sqrt.((x_states[i] .- xmeas).^2 .+ (y_states[i] .- ymeas).^2)./as.localizationDist)

    #vv = Vector(1:convert(Float64, size(densEnsemble,1)))

    #println(size(densEnsemble))
    #println(size(M))
    #println(size(y))
    #println(size(R))
    #println(size(vv))
    #println(size(M*densEnsemble))
    #println(typeof(densEnsemble))
    #println(typeof(M))
    #println(typeof(Rvec))
    #println(typeof(y))
    #println(typeof(vv))
    #println(typeof(selectObs(21)))
    #testloc = selectObs(210)
    #println(string(minimum(testloc))*"   "*string(maximum(testloc)))

    #X_upd, X_upd_mean = DataAssim.local_ESTKF(densEnsemble, M, y, Rvec, vv, selectObs)
    #meanField2 = reshape(X_upd_mean, dimensions[1], dimensions[2])
    X_upd = enKF(densEnsemble, M, xloc, y, Rval, as) # Get corrected ensemble matrix X_upd
    
    #println("densEnsemble: "*string(size(densEnsemble)))
    #println("X_upd: "*string(size(X_upd)))

    X_upd = max.(0.0, real(X_upd))

    #A = densEnsemble - (1/N)*densEnsemble*ones(N,1)*ones(1,N)
    #stdA = std(A, dims=2)
    #display(plot(stdA))

    # Show before and after in density plane:
    twinField = reshape(densTwin[1:nPos], dimensions[1], dimensions[2])
    meanField = reshape(mean(densEnsemble[1:nPos,:], dims=2), dimensions[1], dimensions[2])
    updFieldPre = reshape(mean(X_upd[1:nPos,:], dims=2), dimensions[1], dimensions[2])
    devi = updFieldPre-meanField
    #display(plot(heatmap(meanField), heatmap(twinField), heatmap(updField), heatmap(devi)))
    
    
    # Apply corrections for density for each ensemble member:
    if as.resampleAll # Full resampling strategy

        #aSums = zeros(as.N, 2)   
        for i=1:as.N

            # Compute the deviation field for this member:
            origField = reshape(densEnsemble[1:nPos,i], dimensions[1], dimensions[2])
            densityField = reshape(X_upd[1:nPos,i], dimensions[1], dimensions[2])
            devi = densityField - origField

            # Get the IBM for this ensemble member:
            indsArray = ensemble[i]
            updArray = indsArray
            maxpasses = 100
            doWrite = i==1

            energyField, xrng, yrng = computeAverageEnergyField(updArray, xlim, ylim, dxy, 0.0)
            
            updArray = applyCorrectionsResample(copy(updArray), densityField, energyField, xlim, ylim, dxy, doWrite)            
            ensemble[i] = updArray
        end

    else # Adjustment strategy

        # Call Sinkhorn OT algorithm for all ensemble members:
        ot = applyCorrectionsSinkhornParallel(X_upd[1:nPos,:], densEnsemble[1:nPos,:], dimensions, dxy, true)

        # Then apply results for each ensemble member:
        for i=1:as.N
            #println("i="*string(i))
            # Compute the deviation field for this member:
            origField = reshape(densEnsemble[1:nPos,i], dimensions[1], dimensions[2])
            densityField = reshape(X_upd[1:nPos,i], dimensions[1], dimensions[2])
            
            # Get the IBM for this ensemble member:
            indsArray = ensemble[i]
            updArray = indsArray
        
            doWrite = i==1
            
            # Apply Sinkhorn individually for ensemble members (slow!):
            #updArray, stats, origField = applyCorrectionsSinkhorn(copy(updArray), densityField, origField, xlim, ylim, dxy, doWrite)
            
            # Move individuals according to transport plan:
            #println("Updating member "*string(i)*" based on sinkhorn transport.")
            updArray, origField = applyCorrectionsFromSinkhorn(copy(indsArray), ot[:,:,i], origField, dimensions, xlim, ylim, dxy, doWrite)
            
            # Make size adjustments to match target densities:
            #println("Adjusting sizes")
            updArray, origField, updatedCells = applyCorrectionsResize(copy(updArray), densityField, origField, xlim, ylim, dxy, false)
            ensemble[i] = updArray
            
                            
        end
                
    end
    
    # If activated, apply corrections for speed:
    if as.speedsInStateVec

        for ensI=1:as.N
            # Get corrected speed field for this ensemble member:
            fieldU = reshape(X_upd[stateOffsetSpeed*nPos+1:(stateOffsetSpeed+1)*nPos,ensI], dimensions[1], dimensions[2])
            fieldV = reshape(X_upd[(stateOffsetSpeed+1)*nPos+1:(stateOffsetSpeed+2)*nPos,ensI], dimensions[1], dimensions[2])

            # Iterate over all individuals, find cell and correct speed:
            indsArray = ensemble[ensI]
            for i=1:length(indsArray)
                ix = Int(floor((indsArray[i].x - xlim[1])/dxy))
                iy = Int(floor((indsArray[i].y - ylim[1])/dxy))
                if ix>=1 && iy>=1 && ix<=size(fieldU,1) && iy<=size(fieldU,2)
                    indsArray[i].v_x = fieldU[ix,iy]
                    indsArray[i].v_y = fieldV[ix,iy]
                end
            end
        end

    end

    # If activated, apply corrections for the food field:
    if as.foodInStateVec
        for ensI=1:as.N
            X_fld[:,:,ensI] = reshape(X_upd[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos,ensI], dimensions[1], dimensions[2])
        end
    end

    if doPlot
        updDensEnsemble = zeros(Float64,length(xrng)*length(yrng), as.N)
        for ensI = 1:as.N
            indsArray = ensemble[ensI]
            densityField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)
            updDensEnsemble[:,ensI] = reshape(densityField, length(densityField), 1)
        end    
        updFieldPost = reshape(mean(updDensEnsemble, dims=2), dimensions[1], dimensions[2])

        ind = getStateIndex(dimensions,15, 20)
        covMat = getCorrelationMatrix(dimensions, densEnsemble, ind)

        
        #heatmap(reshape(xloc[:,1],dimensions[1], dimensions[2]),title="Localization matrix"), 

        display(plot(heatmap(twinField, title="Twin field"), heatmap(meanField,title="Orig mean"), 
        heatmap(twinField-meanField,title="Twin - orig"), 
        heatmap(covMat,title="Pre covariance"), 
        heatmap(updFieldPre,title="Updated pre-IBM"), heatmap(updFieldPre-twinField,title="Twin - pre-IBM"),
        heatmap(updFieldPost,title="Updated mean"), heatmap(updFieldPre-updFieldPost,title="Pre-IBM - post-IBM", clim=(-5, 5)), 
        heatmap(updFieldPost-meanField, title="Updated - orig")))

        if as.speedsInStateVec
            # Store ens member 1 u after:
            tmpU, tmpV = computeAverageSpeedField(ensemble[1], xlim, ylim, dxy, 0.0)
            writedlm("C:/temp/Uafter.csv", tmpU, ',')
            writedlm("C:/temp/Vafter.csv", tmpV, ',')
        end

        if as.foodInStateVec
            # Store ens member 1 food after:
            writedlm("C:/temp/Xafter.csv", X_fld[:,:,1], ',')
        end
    end
    
    # Return the updated ensemble and the corrected density field straight from the EnKF:
    return ensemble, X_fld, updFieldPre
end
