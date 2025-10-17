
function ibmAssimilation(as, ensemble, xlim, ylim, dxy, doPlot)

    Ndim = as.N + 1

    # Define ranges:
    xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
    yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
    dimensions = [length(xrng), length(yrng)]
    nPos = length(xrng)*length(yrng)
    statesPerPos = 1
    if as.speedsInStateVec
        statesPerPos = 3
    end
    densEnsemble = zeros(Float64,statesPerPos*nPos, as.N)
    densTwin = zeros(Float64, statesPerPos*nPos, 1)
    #println("defined densEnsemble: ", size(densEnsemble))

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
                densEnsemble[nPos+1:2*nPos,ensI] = reshape(tmpU, length(densityField), 1)
                densEnsemble[2*nPos+1:3*nPos,ensI] = reshape(tmpV, length(densityField), 1)
            else
                densTwin[nPos+1:2*nPos] = reshape(tmpU, length(densityField), 1)
                densTwin[2*nPos+1:3*nPos] = reshape(tmpV, length(densityField), 1)
            end

        end
    end    

    println("Member 1 sum abs deviation: ", sum(abs.(densEnsemble[:,1]-densTwin)) )

    # Ensemble Kalman filter:
    M = getM(dimensions, as) # Measurement model
    xloc = getLocMatrix(dimensions, M, 4)
    y = M*densTwin # Measurement vector based on twin
    Rval = 2.0 # Assumed measurement uncertainty
    X_upd = enKF(densEnsemble, M, xloc, y, Rval) # Get corrected ensemble matrix X_upd

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

        for i=1:as.N

            # Compute the deviation field for this member:
            origField = reshape(densEnsemble[:,i], dimensions[1], dimensions[2])
            densityField = reshape(X_upd[:,i], dimensions[1], dimensions[2])
            devi = densityField - origField

            # Get the IBM for this ensemble member:
            indsArray = ensemble[i]
            updArray = indsArray
        
            doWrite = i==1

            #updArray, origField = applyCorrectionsMoveDirectMultiple(copy(updArray), densityField, origField, xlim, ylim, dxy, maxpasses, doWrite)
            #updArray, origField, updatedCells = applyCorrectionsResize(copy(updArray), densityField, origField, xlim, ylim, dxy, false)
            #if i==1
                #println("Warning, calling only sinkhorn for ensemble member 1!!!")

            updArray, origField = applyCorrectionsSinkhorn(copy(updArray), densityField, origField, xlim, ylim, dxy, doWrite)
            updArray, origField, updatedCells = applyCorrectionsResize(copy(updArray), densityField, origField, xlim, ylim, dxy, false)
            ensemble[i] = updArray
            #end
    
        end
        
    end
    
    # If activatd, apply corrections for speed:
    if as.speedsInStateVec

        for ensI=1:as.N
            # Get corrected speed field for this ensemble member:
            fieldU = reshape(X_upd[nPos+1:2*nPos,ensI], dimensions[1], dimensions[2])
            fieldV = reshape(X_upd[2*nPos+1:3*nPos,ensI], dimensions[1], dimensions[2])

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


    if doPlot
        updDensEnsemble = zeros(Float64,length(xrng)*length(yrng), as.N)
        for ensI = 1:as.N
            indsArray = ensemble[ensI]
            densityField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)
            updDensEnsemble[:,ensI] = reshape(densityField, length(densityField), 1)
        end    
        updFieldPost = reshape(mean(updDensEnsemble, dims=2), dimensions[1], dimensions[2])

        ind = getStateIndex(dimensions,20, 20)
        covMat = getCorrelationMatrix(dimensions, densEnsemble, ind)
        
        display(plot(heatmap(twinField, title="Twin field"), heatmap(meanField,title="Orig mean"), heatmap(twinField-meanField,title="Twin - orig"), 
        #heatmap(reshape(xloc[:,1],dimensions[1], dimensions[2]),title="Localization matrix"), 
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
    end
    

    return ensemble
end
