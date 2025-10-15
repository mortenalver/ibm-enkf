
function getStateIndex(dimensions, i, j)
    return (j-1)*dimensions[1] + i
end

function getCorrelationMatrix(dimensions, X, k)
    szx = size(X)
    myVec = X[k,:]
    cx = zeros(Float64, dimensions[1], dimensions[2])
    for i=1:dimensions[1]
        for j=1:dimensions[2]
            idx = getStateIndex(dimensions, i, j)
            cx[i,j] = cov(myVec, X[idx,:])
        end
    end
    return cx
end


function getRandomField(dims, intensity, ndots, dotRng)
    field = zeros(dims[1], dims[2])
    dots = zeros(ndots, 3)
    for i=1:ndots
        dots[i,:] = [1.0*dims[1]*rand() 1.0*dims[2]*rand() randn()]
    end
    
    for i=1:dims[1]
        for j=1:dims[2]
            for k=1:ndots
                distVec = [i-dots[k,1] j-dots[k,2]]./dotRng
                dist = sqrt(distVec[1]*distVec[1] + distVec[2]*distVec[2])
                effect = exp(-(dist*dist))
                field[i,j] = field[i,j] + intensity*dots[k,3]*effect
            end
        end
    end

    return field
end