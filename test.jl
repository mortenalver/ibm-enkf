
Xf = randn(Float64, 5, 25)
M = [1 0 0 0 0; 0 1 0 0 0]
y = [0.5; 0.5]
R = [1, 1]
vv = Vector(1:convert(Float64, 5))
selectObs(i) = Vector(ones(Float64, length(y)))

println(selectObs(1))
println(size(selectObs(1)))

println(typeof(Xf))
println(typeof(M))
println(typeof(R))
println(typeof(y))
println(typeof(vv))
println(typeof(selectObs(21)))

X_upd, X_upd_mean = DataAssim.local_ESTKF(Xf, M, y, R, vv, selectObs)

println(size(X_upd))
    