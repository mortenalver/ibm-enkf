
function compute_gradient(matrix::Array{Float64, 2})
    # ChatGPT generated function that computes internal gradients of matrix
    rows, cols = size(matrix)
    grad_x = zeros(Float64, rows, cols)
    grad_y = zeros(Float64, rows, cols)

    for i in 2:rows-1
        for j in 2:cols-1
            grad_x[i, j] = (matrix[i+1, j] - matrix[i-1, j]) / 2.0
            grad_y[i, j] = (matrix[i, j+1] - matrix[i, j-1]) / 2.0
        end
    end

    return grad_x, grad_y
end
    
function getIndex(x, y, xlim, ylim, dxy)
    return Int(floor((x - xlim[1])/dxy)), Int(floor((y - ylim[1])/dxy))
end

function applyCorrectionsResize(indsArray, correctedField, origField, xlim, ylim, dxy, doPrint)
    aSumBefore = sum(abs.(correctedField-origField))
    updatedField = origField
    densCorr = correctedField-updatedField
    updatedCells = zeros(size(densCorr))
    adjInd = 0
    testX = 10
    testY = 12
    printedCorr = false
    for ind in indsArray
        # Find which cell this individual is in:
        ix, iy = getIndex(ind.x, ind.y, xlim, ylim, dxy)
        if ix>=1 && iy>=1 && ix<=size(origField,1) && iy<=size(origField,2)
            if origField[ix, iy] > 0
                # Check the relative adjustment needed here:
                relDensCorr = correctedField[ix, iy]/origField[ix, iy]
                # Set limits on the correction: positive and less than a certain maximum value
                relDensCorr = min(5.0, max(0.01, relDensCorr))

                if doPrint && !printedCorr && ix==testX && iy==testY
                    println("Corr fac: ", relDensCorr)
                    printedCorr = true
                end
                # Adjust individual's n value:
                ind.n = ind.n*relDensCorr        
                adjInd = adjInd+1
                #if relDensCorr < 1e-3
                #    println("Warning: zeroing out individual")
                #end
                updatedCells[ix, iy] = updatedCells[ix, iy]+1 # Mark that we have operated in this cell
            end
        end
    end
    updatedField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)

    if doPrint
        println("Adjusted/Total inds: ", adjInd, " / ", length(indsArray))
    end
    
    return indsArray, updatedField, updatedCells
end

function applyCorrectionsMove(indsArray, correctedField, origField, xlim, ylim, dxy)
    
    aSumBefore = sum(abs.(correctedField-origField))
    
    updatedField = origField
    for passes = 1:1
        prevSum = sum(abs.(correctedField-updatedField))
        densCorr = correctedField-updatedField
        gradX, gradY = compute_gradient(densCorr)
        for ind in indsArray
            # Find which cell this individual is in:
            ix, iy = getIndex(ind.x, ind.y, xlim, ylim, dxy)
            # Check if we are to adjust the density of this cell downwards:
            if densCorr[ix, iy] < 0
                # If so, move the individual in the direction of its gradient:
                grad = [gradX[ix, iy], gradY[ix, iy]]
                gradLen = sqrt(grad[1]*grad[1] + grad[2]*grad[2])
                maxGradLen = 10
                if gradLen>maxGradLen
                    grad = grad * maxGradLen/gradLen
                end
                if gradLen > 0
                    ind.x = ind.x + 0.1*grad[1]/gradLen
                    ind.y = ind.y + 0.1*grad[2]/gradLen
                end
            end
        end
        updatedField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)
        newSum = sum(abs.(correctedField-updatedField))    
    end
    aSumAfter = sum(abs.(correctedField-updatedField))
    #println("Abs sum change: ", aSumAfter-aSumBefore)
    return indsArray, aSumAfter-aSumBefore
end


function applyCorrectionsMoveDirectMultiple(updArray, densityField, origField, xlim, ylim, dxy, maxpasses, doWrite)
    goOn = true
    count = 0
    aSumAfter = 0
    while goOn
        testArray, aSumAfter, nmoves, updatedCells = applyCorrectionsMoveDirect(copy(updArray), densityField, origField, xlim, ylim, dxy)
        
        if doWrite
            println("Count ", count, ", aSumAfter = ", aSumAfter, ", nmoves=", nmoves)
        end
        if count > maxpasses || nmoves==0 || aSumAfter > 0
            goOn = false
        else
            updArray = testArray
            origField, xrng, yrng = computeDensityField(updArray, xlim, ylim, dxy)
            count = count+1
        end

        # if !goOn
        #     if doPlot && i==1
        #         updDens = zeros(Float64,length(xrng)*length(yrng))
        #         updDensityField, xrng, yrng = computeDensityField(updArray, xlim, ylim, dxy)
        #         println("Devi after at (10,12): ", densityField[10,12]/updDensityField[10,12])
        #         display(plot(heatmap(updDensityField, title="Adjusted field"), 
        #             heatmap(densityField-updDensityField,title="Pre-IBM - post-IBM",clim=(-1,1)), 
        #             heatmap(updatedCells, title="Updated cells", clim=(0,1))))
        #     end
        # end

    end

    return updArray, origField
end

function applyCorrectionsMoveDirect(indsArray, correctedField, origField, xlim, ylim, dxy)
    nmoves = 0
    #aSumBefore = sum(abs.(correctedField-origField))
    aSumBefore = sum((correctedField-origField).^2)

    updatedCells = zeros(size(origField))

    updatedField = origField

    densCorr = correctedField-updatedField
    for ind in indsArray
        # Find which cell this individual is in:
        ix, iy = getIndex(ind.x, ind.y, xlim, ylim, dxy)
        # Check if we are to adjust the density of this cell downwards:
        if ix>=1 && iy>=1 && ix<=size(densCorr,1) && iy<=size(densCorr,2)
            if densCorr[ix, iy] < -ind.n #0
                # If so, look for a nearby cell that needs more individuals:
                #nearOffs = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1;]
                nearOffs = [1 0; 0 1; -1 0; 0 -1; 1 1; 1 -1; -1 1; -1 -1;
                    -2 -1; -2 0; -2 1; 2 -1; 2 0; 2 1; -1 2; 0 2; 1 2; -1 -2; 0 -2; 1 -2;
                    -2 -2; -2 2; 2 -2; 2 2;
                    -3 -1; -3 0; -3 1; 3 -1; 3 0; 3 1; -1 -3; 0 -3; 1 -3; -1 3; 0 3; 1 3]
                for ii in shuffle(1:size(nearOffs,1))
                    coord = [ix;iy] + nearOffs[ii,:]
                    if coord[1]>=1 && coord[2]>=1 && coord[1]<=size(densCorr,1) && coord[2]<=size(densCorr,2)
                        if densCorr[coord[1], coord[2]] >= ind.n
                            if updatedField[coord[1], coord[2]] == 0 # Only move to empty cell
                                # Move to that cell:
                                ind.x = ind.x + nearOffs[ii,1]
                                ind.y = ind.y + nearOffs[ii,2]
                                densCorr[ix,iy] = densCorr[ix,iy] + ind.n
                                densCorr[coord[1],coord[2]] = densCorr[coord[1],coord[2]] - ind.n
                                nmoves = nmoves+1
                                updatedCells[ix, iy] = 1
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    updatedField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)


    #aSumAfter = sum(abs.(correctedField-updatedField))
    aSumAfter = sum((correctedField-updatedField).^2)
    
    return indsArray, aSumAfter-aSumBefore, nmoves, updatedCells
end
