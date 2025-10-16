
function applyCorrectionsSinkhorn(indsArray, correctedField, origField, xlim, ylim, dxy, doWrite)
    #println("Start apply sinkhorn")
    # Subsample fields to reduce complexity of optimization problem:
    subs = 2
    #println(subs)
    corr = zeros(Float64, convert(Int, ceil(size(origField,1)/subs)), convert(Int, ceil(size(origField,2)/subs)))
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
            C[i,j] = dist.^0.8
        end
    end
    if doWrite
        println("Calling sinkhorn algorithm. sizes: "*string(size(corr))*" "*string(size(a))*" "*string(size(C)))
    end
    ot = sinkhorn(a, b, C, 0.95, maxiter=100000)#, alg=SinkhornGibbs())
    if ~isnan(ot[1,1])
        #println("Moving inds....")
        # Go through each individual, find its cell and the transport distribution for that cell:
        for i=1:length(indsArray)
            #if doWrite
            # Find cell in subsampled array:
            ix = Int(floor((indsArray[i].x - xlim[1])/(subs*dxy)))
            iy = Int(floor((indsArray[i].y - ylim[1])/(subs*dxy)))
            if ix>=1 && iy>=1 && ix<=size(orig,1) && iy<=size(orig,2)
                indx = size(orig,1)*(iy-1)+ix
                #println("i="*string(i)*" x="*string(indsArray[i].x)*" y="*string(indsArray[i].y)*" ix="*string(ix)*" iy="*string(iy)*" indx="*string(indx))
                toVals = ot[indx,:]
                if sum(toVals)==0.0
                    continue
                end
                permu = sortperm(toVals, rev=true)
                
                # Find relative probabilities of the first 8 target cells:
                relProb = toVals[permu[1:8]] ./ sum(toVals[permu[1:8]])
                cumProb = cumsum(relProb)
                #println("First to cells: "*string(permu[1:8]))
                #println("relProb: "*string(relProb))
                #println("cumProb: "*string(cumProb))

                # Choose one of those cells randomly according to their probabilities:
                cellNo = searchsortedfirst(cumProb, rand())
                #println("cellNo: "*string(cellNo))
                
                # Find cell location of the cell to move to:
                toInd = permu[cellNo]
                if toInd!=indx
                    cellI = cind[toInd][1]
                    cellJ = cind[toInd][2]
                    # if abs(cellI-ix) > 2 || abs(cellJ-iy) > 2
                    #     println("Move to: "*string(cellI)*" , "*string(cellJ)*": deltaX="*string(subs*dxy*(cellI-ix))*": deltaY="*string(subs*dxy*(cellJ-iy))*" sumToVals="*string(sum(toVals)))
                    #     println(string(cellNo)*" "*string(permu[1:8])*" "*string(relProb))
                    # end
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

    updatedField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)

    if doWrite
        println("Writing result to file")
        writedlm("C:/temp/orig.csv", orig, ',')
        writedlm("C:/temp/orig0.csv", origField, ',')
        writedlm("C:/temp/corr.csv", corr, ',')
        writedlm("C:/temp/corr0.csv", correctedField, ',')
        writedlm("C:/temp/OT.csv", ot, ',')
        writedlm("C:/temp/postdens.csv", updatedField,',')
    end


    return indsArray, updatedField
end