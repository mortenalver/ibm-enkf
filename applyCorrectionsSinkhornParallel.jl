
# Function to calculate RMS for a vector of numbers
function compute_rms(data)
    @assert length(data) > 0 "Input array cannot be empty"
    return sqrt(sum((data .- mean(data)).^2) / length(data))
end

function applyCorrectionsSinkhornParallel(correctedField, origField, dims, dxy, doWrite)
    println("Start apply sinkhorn")
    
    a = origField#[:,1:10]
    b = copy(correctedField)#[:,1:10])
    for i=1:size(a,2)
        b[:,i] = b[:,i] ./ (sum(b[:,i])/sum(a[:,i])) # Normalize so their sums are equal
    end

    # Set up cost matrix:
    # if doWrite
    #     println("Setting up cost matrix")
    # end
    C = zeros(size(a,1), size(a,1))
    cind = CartesianIndices((dims[1], dims[2]))
    for i=1:size(a,1)
        i1 = cind[i][1]
        j1 = cind[i][2]
        for j=1:size(a,1)
            i2 = cind[j][1]
            j2 = cind[j][2]

            distVec = [i1-i2 j1-j2]
            dist = distVec[1]*distVec[1]+distVec[2]*distVec[2]
            C[i,j] = dist#.^0.8
        end
    end
    if doWrite
        println("Calling sinkhorn algorithm. sizes: "*string(size(a))*" "*string(size(b))*" "*string(size(C)))
    end

    aa = copy(a)
    println(size(aa))
    # Call sinkhorn algotithm:
    eps = 1.9
    ot = sinkhorn(aa, b, C, eps, 
        SinkhornEpsilonScaling(
            SinkhornGibbs();
            factor=0.55,#1//2,  # 0.65: 595 sek
            steps=2,
        );
        atol=1e-6,
        maxiter=50_000)

    # println(size(ot))
    # for i=1:size(ot,3)
    #     # Compute row sums to check against source distribution:
    #     rSums = sum(ot[:,:,i],dims=2)
    #     rowRMS = compute_rms(rSums-a[:,i])
    #     cSums = sum(ot[:,:,i],dims=1)'
    #     colRMS = compute_rms(cSums-b[:,i])
    #     println(string(i)*": Row RMS: "*string(rowRMS)*", Col RMS: "*string(colRMS))
    # end

    #optCost = sinkhorn2(a, b, C, 1.0, maxiter=200000, plan=ot)
    
    return ot
end


function applyCorrectionsFromSinkhorn(indsArray, ot, orig, dims, xlim, ylim, dxy, doWrite)

    if ~isnan(ot[1,1])
        cind = CartesianIndices((dims[1], dims[2]))
        # Go through each individual, find its cell and the transport distribution for that cell:
        for i=1:length(indsArray)
            #if doWrite
            # Find cell in array:
            ix = Int(floor((indsArray[i].x - xlim[1])/dxy))
            iy = Int(floor((indsArray[i].y - ylim[1])/dxy))
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
                    indsArray[i].x += dxy*(cellI-ix)
                    indsArray[i].y += dxy*(cellJ-iy)
                end
            end

        end
        #println("...done")
    else
        println("Sinkhorn algorithm returned NaN")
    end

    updatedField, xrng, yrng = computeDensityField(indsArray, xlim, ylim, dxy)


    return indsArray, updatedField

end