
function ibmAssimilation(as, ensemble, xlim, ylim, dxy, doPlot)

    Ndim = as.N + 1

    # Define ranges:
    xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
    yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
    dimensions = [length(xrng), length(yrng)]

    densEnsemble = zeros(Float64,length(xrng)*length(yrng), as.N)
    densTwin = zeros(Float64, length(xrng)*length(yrng), 1)
    #println("defined densEnsemble: ", size(densEnsemble))

    for ensI = 1:Ndim
        indsArray = ensemble[ensI]
        densityField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)

        if ensI <= as.N
            densEnsemble[:,ensI] = reshape(densityField, length(densityField), 1)
        else
            densTwin = reshape(densityField, length(densityField), 1)
        end
    end    

    println("Member 1 sum abs deviation: ", sum(abs.(densEnsemble[:,1]-densTwin)) )

    # Ensemble Kalman filter:
    M = getM(dimensions) # Measurement model
    xloc = getLocMatrix(dimensions, M, 4)
    y = M*densTwin # Measurement vector based on twin
    Rval = 2.0 # Assumed measurement uncertainty
    X_upd = enKF(densEnsemble, M, xloc, y, Rval) # Get corrected ensemble matrix X_upd

    #A = densEnsemble - (1/N)*densEnsemble*ones(N,1)*ones(1,N)
    #stdA = std(A, dims=2)
    #display(plot(stdA))

    # Show before and after in density plane:
    twinField = reshape(densTwin, dimensions[1], dimensions[2])
    meanField = reshape(mean(densEnsemble, dims=2), dimensions[1], dimensions[2])
    updFieldPre = reshape(mean(X_upd, dims=2), dimensions[1], dimensions[2])
    devi = updFieldPre-meanField
    #display(plot(heatmap(meanField), heatmap(twinField), heatmap(updField), heatmap(devi)))
    
    
    # Apply corrections for each ensemble member:
    aSums = zeros(as.N, 2)   
    for i=1:as.N

        # Compute the deviation field for this member:
        origField = reshape(densEnsemble[:,i], dimensions[1], dimensions[2])
        densityField = reshape(X_upd[:,i], dimensions[1], dimensions[2])
        devi = densityField - origField

        # Get the IBM for this ensemble member:
        indsArray = ensemble[i]
        updArray = indsArray
        maxpasses = 100
        doWrite = i==1

        if as.resampleAll # Full resampling strategy

            energyField, xrng, yrng = computeAverageEnergyField(updArray, xlim, ylim, dxy, 0.0)
            
            updArray = applyCorrectionsResample(copy(updArray), densityField, energyField, xlim, ylim, dxy, doWrite)            
            ensemble[i] = updArray

        else # Adjustment strategy
            #updArray, origField = applyCorrectionsMoveDirectMultiple(copy(updArray), densityField, origField, xlim, ylim, dxy, maxpasses, doWrite)
            #updArray, origField, updatedCells = applyCorrectionsResize(copy(updArray), densityField, origField, xlim, ylim, dxy, false)
            #if i==1
                #println("Warning, calling only sinkhorn for ensemble member 1!!!")

                updArray, origField = applyCorrectionsSinkhorn(copy(updArray), densityField, origField, xlim, ylim, dxy, doWrite)
                updArray, origField, updatedCells = applyCorrectionsResize(copy(updArray), densityField, origField, xlim, ylim, dxy, false)
                ensemble[i] = updArray
            #end
                    
        end

        #ensemble[i] = updArray
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
    end
    

    return ensemble
end
