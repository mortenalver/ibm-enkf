
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


function applyCorrectionsSinkhorn(indsArray, correctedField, origField, xlim, ylim, dxy, doWrite)
    #println("Start apply sinkhorn")
    # Subsample fields to reduce complexity of optimization problem:
    subs = 2
    corr = zeros(Float64, convert(Int, ceil(size(origField,1)/2)), convert(Int, ceil(size(origField,2)/2)))
    orig = zeros(Float64, size(corr))
    for i=1:size(orig,1)
        for j=1:size(orig,2)
            orig[i,j] = sum(origField[(subs*i-1):min(size(origField,1), subs*i), (subs*j-1):min(size(origField,2), subs*j)])
            corr[i,j] = sum(correctedField[(subs*i-1):min(size(origField,1), subs*i), (subs*j-1):min(size(origField,2), subs*j)])
        end
    end
    #println(size(corr))

    # Flatten original and corrected fields:
    a = reshape(orig, length(orig))
    b = reshape(corr, length(corr))
    b = b ./ (sum(b)/sum(a)) # Normalize so their sums are equal
    # Set up cost matrix:
    if doWrite
        println("Setting up cost matrix")
    end
    C = zeros(length(a), length(a))
    cind = CartesianIndices(size(orig))
    for i=1:length(a)
        i1 = cind[i][1]
        j1 = cind[i][2]
        for j=1:length(a)
            i2 = cind[j][1]
            j2 = cind[j][2]

            distVec = [i1-i2 j1-j2]
            dist = distVec[1]*distVec[1]+distVec[2]*distVec[2]
            C[i,j] = dist
        end
    end
    if doWrite
        println("Calling sinkhorn algorithm. sizes: "*string(size(corr))*" "*string(size(a))*" "*string(size(C)))
    end
    ot = sinkhorn(a, b, C, 0.9, maxiter=20000)#, alg=SinkhornGibbs())
    if ~isnan(ot[1,1])
        #println("Moving inds....")
        # Go through each individual, find its cell and the transport distribution for that cell:
        for i=1:length(indsArray)
            #if doWrite
            # Find cell in subsampled array:
            ix = Int(floor((indsArray[i].x - xlim[1])/(2.0*dxy)))
            iy = Int(floor((indsArray[i].y - ylim[1])/(2.0*dxy)))
            if ix>=1 && iy>=1 && ix<=size(orig,1) && iy<=size(orig,2)
                indx = size(orig,1)*(iy-1)+ix
                #println("i="*string(i)*" x="*string(indsArray[i].x)*" y="*string(indsArray[i].y)*" ix="*string(ix)*" iy="*string(iy)*" indx="*string(indx))
                toVals = ot[indx,:]
                permu = sortperm(toVals, rev=true)
                
                # Find relative probabilities of the first 5 target cells:
                relProb = toVals[permu[1:5]] ./ sum(toVals[permu[1:5]])
                cumProb = cumsum(relProb)
                #println("First to cells: "*string(permu[1:5]))
                #println("relProb: "*string(relProb))
                #println("cumProb: "*string(cumProb))
                cellNo = searchsortedfirst(cumProb, rand())
                #println("cellNo: "*string(cellNo))
                
                # Find cell location of the cell to move to:
                toInd = permu[cellNo]
                if toInd!=indx
                    cellI = cind[toInd][1]
                    cellJ = cind[toInd][2]
                    #println("Move to: "*string(cellI)*" , "*string(cellJ)*": deltaX="*string(subs*dxy*(cellI-ix))*": deltaY="*string(subs*dxy*(cellJ-iy)))
                    # Update position:
                    indsArray[i].x += subs*dxy*(cellI-ix)
                    indsArray[i].y += subs*dxy*(cellJ-iy)
                end
            end

        end
        #println("...done")
    else
        println("Sinkhorn algorithm returned NaN")
    end

    if doWrite
        println("Writing result to file")
        writedlm("C:/temp/orig.csv", orig, ',')
        writedlm("C:/temp/orig0.csv", origField, ',')
        writedlm("C:/temp/corr.csv", corr, ',')
        writedlm("C:/temp/corr0.csv", correctedField, ',')
        writedlm("C:/temp/OT.csv", ot, ',')
    end

    updatedField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)

    return indsArray, updatedField
end