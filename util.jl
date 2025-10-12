
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