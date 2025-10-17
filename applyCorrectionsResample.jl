
function applyCorrectionsResample(indsArray, correctedField, energyField, xlim, ylim, dxy, doPrint)
    # Adapt individual model by resampling all individuals randomly.

    # Make a map of individuals vs. grid indexes:
    totInd = 0.0
    idx = 0
    indMap = zeros(Float64, length(correctedField))
    #for i=1:length(indMap)
    for j=1:size(correctedField,2)
        for i=1:size(correctedField,1)        
            idx += 1
            totInd = totInd + correctedField[i,j]
            indMap[idx] = totInd
        end
    end

    # Then roll uniform random numbers once per individual,
    # choose cell and move the individual:
    cind = CartesianIndices(size(correctedField))
    for i in eachindex(indsArray)
        ind = indsArray[i]
        randNum = totInd*rand()
        cellNo = searchsortedfirst(indMap, randNum)
        cellI = cind[cellNo][1]
        cellJ = cind[cellNo][2]
        #xold = ind.x
        ind.x = (cellI+rand())*dxy
        ind.y = (cellJ+rand())*dxy

        # Update energy of this individual to the mean energy of the cell it moves to:
        ind.E = energyField[cellI,cellJ]

    end

    return indsArray
end
