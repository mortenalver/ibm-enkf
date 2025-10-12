
function enKF(X, M, xloc, d, Rval)
    N = size(X)[2]
    D = zeros(Float64, length(d), N)
    for i=1:N
        D[:,i] = d
    end
    #println("D: ", size(D))
    Rvec = Rval*ones(Float64, length(d))
    R = Diagonal(Rvec)
    #println("R: ", size(R))
    
    #println("N=", N)
    # Compute intermediate matrices:
    MX = M*X
    MA = MX - (1/N)*(MX*ones(N,1))*ones(1,N)
    P = (1/(N-1))*MA*MA' + R
    #println("P: ", size(P))
    # P_loc = (1/(N-1))*((MA*MA').*XlocM) + C_ee;
    A = X - (1/N)*X*ones(N,1)*ones(1,N)
    #println("A: ", size(A))
    K = xloc .* ((1/(1+N))*(A*MA')*inv(P))
    #println("K: ", size(K))
    # K_loc = (1/(1+N))*((A*MA').*Xloc)*inv(P_loc);
    # Compute analysis:
    correction = K*(D - MX)
    #println("corr: ", size(correction))

    return X + correction
end