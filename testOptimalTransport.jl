#import Pkg; Pkg.add("OptimalTransport")
using OptimalTransport 

function testOT()
    #A = [1 0 0 0; 0 0 0 0; 0 0 0 0]
    #B = [0 0 0 0; 0 1.0 1.0 0; 0 0 0 0]
    A = randn(40, 30)
    B = randn(40, 30)
    #A = rand(15, 10)
    #B = rand(15, 10)
    
    B = B ./ (sum(B)/sum(A))
    a = reshape(A, length(A))
    b = reshape(B, length(B))

    C = zeros(length(A), length(A))
    cind = CartesianIndices(size(A))
    for i=1:length(A)
        i1 = cind[i][1]
        j1 = cind[i][2]
        for j=1:length(A)
            i2 = cind[j][1]
            j2 = cind[j][2]

            distVec = [i1-i2 j1-j2]
            dist = distVec[1]*distVec[1]+distVec[2]*distVec[2]
            C[i,j] = dist
        end
    end
    #print(C)
    res = sinkhorn(a, b, C, 0.9, maxiter=100000)#, alg=SinkhornGibbs())

    #res = sinkhorn(
    #    a, b, C, 0.01, alg=SinkhornGibbs();
    #    atol=0, rtol=0.01, check_convergence=10, maxiter=10000,
    #)
    # atol=0, rtol=atol > 0 ? 0 : √eps, check_convergence=10, maxiter=1000,

    #println(res)
    writedlm("C:/temp/OT.csv", res, ',')
end