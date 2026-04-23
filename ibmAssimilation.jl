
function ibmAssimilation(as, ensemble, X_fld, xlim, ylim, dxy, doPlot, t)
    # Handle the assimilation process from derivation of state values
    # via calling the Ensemble Kalman filter to updating the IBM.

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
    #println("statesPerPos: "*string(statesPerPos))
    
    densEnsemble = zeros(Float64,statesPerPos*nPos, as.N)
    densTwin = zeros(Float64, statesPerPos*nPos, 1)

    #println("size densTwin: "*string(size(densTwin)))

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

        if doPlot && ensI==1
            writedlm(storageDir*"dens_before.csv", densityField, ',')
        end

        if ensI <= as.N
            densEnsemble[1:nPos,ensI] = reshape(densityField, length(densityField), 1)
        else
            densTwin[1:nPos] = reshape(densityField, length(densityField), 1)
        end

        # Add speeds to state vector if we are supposed to:
        if as.speedsInStateVec
            tmpU, tmpV = computeAverageSpeedField(indsArray, xlim, ylim, dxy, 0.0)
        
            if doPlot && ensI==1
                writedlm(storageDir*"Ubefore.csv", tmpU, ',')
                writedlm(storageDir*"Vbefore.csv", tmpV, ',')
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
            writedlm(storageDir*"Xbefore.csv", X_fld[:,:,1], ',')
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

    if doPlot 
        writedlm(storageDir*"X_f_all.csv", densEnsemble, ',')
    end

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
    writedlm(storageDir*"Xloc.csv", xloc, ',')
    writedlm(storageDir*"M.csv", M, ',')

    y = M*densTwin # Measurement vector based on twin
    Rval = as.measVar # Assumed measurement uncertainty
    
    Rvec = Rval*ones(Float64, size(M,1))
    R = Diagonal(Rvec)
    
    


    #X_upd, X_upd_mean = DataAssim.ESTKF(densEnsemble, M*densEnsemble, y, R, M)

    #selectObs(i) = exp.(-((x_states[i] .- xmeas).^2 .+ (y_states[i] .- ymeas).^2)./(as.localizationDist.^2))
    #vv = Vector(1:convert(Float64, size(densEnsemble,1)))

    #println(size(densEnsemble))
    #println(size(M))
    #println(size(y))
    #println(size(R))
    #println(size(vv))
    #println(size(M*densEnsemble))
    #testloc = selectObs(210)
    #println(string(minimum(testloc))*"   "*string(maximum(testloc)))
    # function $method(Xf,HXf,y,R,H; debug = false, tolerance=1e-10)
    #X_upd, X_upd_mean = DataAssim.ESTKF(densEnsemble, M*densEnsemble, y, R, M)
    #X_upd, X_upd_mean = DataAssim.local_ESTKF(densEnsemble, M, y, Rvec, vv, selectObs)
    #meanField2 = reshape(X_upd_mean, dimensions[1], dimensions[2])

    if !as.anamorphicTransform
        X_upd = enKF(densEnsemble, M, xloc, y, Rval, Rvec, as) # Get corrected ensemble matrix X_upd
    else
        densEnsemble_trf = copy(densEnsemble)
        # Initialize transform:
        trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp = initTransform(vec(densEnsemble_trf[1:nPos,:]), 150)
        println("trf_x: "*string(trf_x))
        println("trf_xl: "*string(trf_xl))
        println("trf_xp: "*string(trf_xp))
        println("trf_yl: "*string(trf_yl))
        println("trf_yp: "*string(trf_yp))

        # Transform density values:
        for i=1:as.N
            trfVal = transform(densEnsemble_trf[1:nPos,i], trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
            densEnsemble_trf[1:nPos,i] = trfVal
        end

        # Transform measurement series (only density values):
        y_trf = y
        if as.measureFood
            n_meas = length(y)
            n_trf = Int(round(n_meas/2))
            y_trf[1:n_trf] = transform(y_trf[1:n_trf], trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
        else
            y_trf[1:n_meas] = transform(y_trf[1:n_meas], trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
        end
        #println(size(Rval))
        for i=1:n_meas
            Rvec[i] = 1.0^2.0 # Reduce value since we are transforming these variables
        end
        X_upd = enKF(densEnsemble_trf, M, xloc, y, Rval, Rvec, as) # Get corrected ensemble matrix X_upd

        X_upd_nontr = copy(X_upd)
        # Iverse transform of density values:
        for i=1:as.N
            trfVal = invTransform(X_upd[1:nPos,i], trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
            X_upd[1:nPos,i] = trfVal
        end

        #trOrig = reshape(densEnsemble_trf[1:nPos,1], dimensions[1], dimensions[2])
        #trUpd = reshape(X_upd_nontr[1:nPos,1], dimensions[1], dimensions[2])
        #invtrUpd = reshape(X_upd[1:nPos,1], dimensions[1], dimensions[2])
        #display(plot(heatmap(trOrig, title="Orig transformed"), heatmap(trUpd,title="Corrected"),
        #    heatmap(invtrUpd,title="Corrected inv-transformed")))
    end
    #println("densEnsemble: "*string(size(densEnsemble)))
    #println("X_upd: "*string(size(X_upd)))

    # Make sure all values are real:
    X_upd = real(X_upd)
    
    allCount = length(X_upd[1:nPos,:])
    szU = size(X_upd[1:nPos,:])
    negCount = 0
    negSum = 0
    for i=1:szU[1]
        for j=1:szU[2]
            if (X_upd[i,j] < 0.0)
                negCount = negCount+1
                negSum = negSum+X_upd[i,j]
            end
        end
    end
    println("Density: Negative / all: "*string(negCount)*" / "*string(allCount)*" , fraction = "*string(float(negCount)/float(allCount))*" , avg value = "*string(negSum/float(negCount)))
    X_upd_part = X_upd[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos,:]
    allCount = length(X_upd_part)
    szU = size(X_upd_part)
    negCount = 0
    negSum = 0
    for i=1:szU[1]
        for j=1:szU[2]
            if (X_upd_part[i,j] < 0.0)
                negCount = negCount+1
                negSum = negSum+X_upd_part[i,j]
            end
        end
    end
    println("Food: Negative / all: "*string(negCount)*" / "*string(allCount)*" , fraction = "*string(float(negCount)/float(allCount))*" , avg value = "*string(negSum/float(negCount)))
    
    # Cut off any negative values for density and food (but not for speeds, which can be negative):
    X_upd[1:nPos,:] = max.(0.0, X_upd[1:nPos,:])
    if as.foodInStateVec
        X_upd[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos,:] = max.(0.0, X_upd[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos,:])
    end

    # Show before and after in density plane:
    twinField = reshape(densTwin[1:nPos], dimensions[1], dimensions[2])
    meanField = reshape(mean(densEnsemble[1:nPos,:], dims=2), dimensions[1], dimensions[2])
    updFieldPre = reshape(mean(X_upd[1:nPos,:], dims=2), dimensions[1], dimensions[2])
    devi = updFieldPre-meanField
    #display(plot(heatmap(meanField), heatmap(twinField), heatmap(updField), heatmap(devi)))
    
    if doPlot
        dddField = reshape(X_upd[1:nPos,1], dimensions[1], dimensions[2])
        writedlm(storageDir*"dens_after.csv", dddField, ',')
    end
    
    # Apply corrections for density for each ensemble member:
    if as.resampleAll # Full resampling strategy

        #aSums = zeros(as.N, 2)   
        Threads.@threads for i=1:as.N

            # Compute the deviation field for this member:
            origField = reshape(densEnsemble[1:nPos,i], dimensions[1], dimensions[2])
            densityField = reshape(X_upd[1:nPos,i], dimensions[1], dimensions[2])
            devi = densityField - origField

            # Get the IBM for this ensemble member:
            indsArray = ensemble[i]
            updArray = indsArray
            
            doWrite = i==1

            energyField, xrng, yrng = computeAverageEnergyField(updArray, xlim, ylim, dxy, 0.0)
            
            updArray = applyCorrectionsResample(copy(updArray), densityField, energyField, xlim, ylim, dxy, doWrite)            
            ensemble[i] = updArray
        end

    else # Adjustment strategy

        # Call Sinkhorn OT algorithm for all ensemble members:
        ot = applyCorrectionsSinkhornParallel(X_upd[1:nPos,:], densEnsemble[1:nPos,:], dimensions, dxy, true)

        # Then apply results for each ensemble member:
        Threads.@threads for i=1:as.N
            #println("i="*string(i))
            # Compute the deviation field for this member:
            origField = reshape(densEnsemble[1:nPos,i], dimensions[1], dimensions[2])
            densityField = reshape(X_upd[1:nPos,i], dimensions[1], dimensions[2])
            
            # Get the IBM for this ensemble member:
            indsArray = ensemble[i]
            updArray = indsArray
        
            doWrite = i==1
            
            # if doWrite
            #     filename = "OT_example_t_"*string(Int(round(t)))*".csv"
            #     println("Writing "*filename)
            #     writedlm(filename, ot[:,:,i], ',')

            # end

            # Apply Sinkhorn individually for ensemble members (slow!):
            #updArray, stats, origField = applyCorrectionsSinkhorn(copy(updArray), densityField, origField, xlim, ylim, dxy, doWrite)
            
            # Move individuals according to transport plan:
            #println("Updating member "*string(i)*" based on sinkhorn transport.")
            updArray, origField = applyCorrectionsFromSinkhorn(copy(indsArray), ot[:,:,i], origField, dimensions, xlim, ylim, dxy,
                                                               as.fuzzySinkhornMoves, doWrite)
            
            # Make size adjustments to match target densities:
            #println("Adjusting sizes")
            updArray, origField, updatedCells = applyCorrectionsResize(copy(updArray), densityField, origField, xlim, ylim, dxy, false)
            ensemble[i] = updArray
            
                            
        end
                
    end
    
    # If activated, apply corrections for speed:
    if as.speedsInStateVec #&& as.resampleAll

        Threads.@threads for ensI=1:as.N
            # Get corrected speed field for this ensemble member:
            fieldU = reshape(X_upd[stateOffsetSpeed*nPos+1:(stateOffsetSpeed+1)*nPos,ensI], dimensions[1], dimensions[2])
            fieldV = reshape(X_upd[(stateOffsetSpeed+1)*nPos+1:(stateOffsetSpeed+2)*nPos,ensI], dimensions[1], dimensions[2])

            # Iterate over all individuals, find cell and correct speed:
            indsArray = ensemble[ensI]
            for i=1:length(indsArray)
                ix = Int(floor((indsArray[i].x - xlim[1])/dxy))
                iy = Int(floor((indsArray[i].y - ylim[1])/dxy))
                if ix>=1 && iy>=1 && ix<=size(fieldU,1) && iy<=size(fieldU,2)
                    #if as.resampleAll
                        indsArray[i].v_x = fieldU[ix,iy]
                        indsArray[i].v_y = fieldV[ix,iy]
                    #else
                    #    upFrac = 0.2
                    #    indsArray[i].v_x = (1.0-upFrac)*indsArray[i].v_x + upFrac*fieldU[ix,iy]
                    #    indsArray[i].v_y = (1.0-upFrac)*indsArray[i].v_y + upFrac*fieldV[ix,iy]
                    #end
                end
            end
        end

    end

    # If activated, apply corrections for the food field:
    if as.foodInStateVec
        Threads.@threads for ensI=1:as.N
            newXfld = reshape(X_upd[stateOffsetFood*nPos+1:(stateOffsetFood+1)*nPos,ensI], dimensions[1], dimensions[2])
            X_fld[:,:,ensI] = min.(max.(newXfld,0.),as.foodLimit)
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
            writedlm(storageDir*"Uafter.csv", tmpU, ',')
            writedlm(storageDir*"Vafter.csv", tmpV, ',')
        end

        if as.foodInStateVec
            # Store ens member 1 food after:
            writedlm(storageDir*"Xafter.csv", X_fld[:,:,1], ',')
        end
    end
    
    # Return the updated ensemble and the corrected density field straight from the EnKF:
    return ensemble, X_fld, updFieldPre
end

