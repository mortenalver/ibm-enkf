
function findnearest(a,x)
    length(a) > 0 || return 0:-1
    r = searchsorted(a,x)
    length(r) > 0 && return r
    last(r) < 1 && return searchsorted(a,a[first(r)])
    first(r) > length(a) && return searchsorted(a,a[last(r)])
    x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
    x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
    return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end

function computeDensityField(indsArray, xlim, ylim, dxy)
    # Compute a field as given by xlim, ylim and dxy that counts all individuals
    # weighted by their N value per grid cell:
    xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
    yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
    #field = zeros(Float64, length(xrng), length(yrng))
    field = zeros(length(xrng), length(yrng))
    for i = 1:length(indsArray)
        #psX = findnearest(xrng, indsArray[i].x)
        #println("1: ", psX[1])
        ix = Int(floor((indsArray[i].x - xlim[1])/dxy))
        iy = Int(floor((indsArray[i].y - ylim[1])/dxy))
        #println("2: ", index)

        #psY = findnearest(yrng, indsArray[i].y)

        if ix>=1 && iy>=1 && ix<=size(field,1) && iy<=size(field,2)
            field[ix, iy] = field[ix, iy] + indsArray[i].n
        end
    end
    return field, xrng, yrng
end

function computeAverageEnergyField(indsArray, xlim, ylim, dxy, eFillVal)
    # Compute a field as given by xlim, ylim and dxy with average E values 
    # of individuals per cell, weighted by the individuals' N values:
    # weighted by their N value per grid cell:
    xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
    yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
    #field = zeros(Float64, length(xrng), length(yrng))
    fieldE = zeros(length(xrng), length(yrng))
    fieldN = zeros(length(xrng), length(yrng))
    for i = 1:length(indsArray)
        ix = Int(floor((indsArray[i].x - xlim[1])/dxy))
        iy = Int(floor((indsArray[i].y - ylim[1])/dxy))
        
        if ix>=1 && iy>=1 && ix<=size(fieldE,1) && iy<=size(fieldE,2)
            # Add n to the density field:
            fieldN[ix, iy] = fieldN[ix, iy] + indsArray[i].n
            # Add n*E to the summed energy field:
            fieldE[ix, iy] = fieldE[ix, iy] + indsArray[i].n * indsArray[i].E
        end
    end
    # Then divide the fummed energy field by the density field, leaving 0 in cells
    # with no individuals:
    fieldE = fieldE ./ fieldN
    for i=1:length(fieldE)
        if isnan(fieldE[i])
            fieldE[i] = eFillVal
        end
    end
    return fieldE, xrng, yrng
end

function computeDensityFieldArrays(x, y, xlim, ylim, dxy)
    #println(xlim)
    xrng = range(start=xlim[1], stop=xlim[2], step=dxy)
    yrng = range(start=ylim[1], stop=ylim[2], step=dxy)
    #field = zeros(Float64, length(xrng), length(yrng))
    field = zeros(length(xrng), length(yrng))
    for i = 1:length(x)
        psX = findnearest(xrng, x[i])
        psY = findnearest(yrng, y[i])
        field[psX[1], psY[1]] = field[psX[1], psY[1]] + 1
    end
    println(size(field))
    return field
end