using Distributions


function initTransform(a, nPercentiles)
    println(typeof(a))
    a_srt = sort(a)
    n_a = length(a_srt)
    trf_x = fill(0.0, nPercentiles)
    trf_y = fill(0.0, nPercentiles)
    for i=1:nPercentiles
        rk = Int(round(float(n_a)*float(i)/float(nPercentiles)))
        trf_x[i] = a_srt[rk]
        trf_y[i] = quantile(Normal(), min(0.99,float(i)/float(nPercentiles)))
    end
    trf_xl = minimum(a)
    trf_xp = maximum(a)
    trf_yl = quantile(Normal(), 1.0/float(nPercentiles)/2.0)
    trf_yp = maximum(trf_y)

    return trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp
end

function transform(a, trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
    p = length(trf_x)
    a_tr = fill(0.0, length(a))
    for i=1:length(a)
        inds = findall(a[i].>trf_x)
        if length(inds)==0
            if a[i] > trf_xl
                a_tr[i] = trf_yl + ((trf_y[1]-trf_yl)/(trf_x[1]-trf_xl))*(a[i]-trf_xl)
            else
                a_tr[i] = trf_yl
            end
        elseif maximum(inds)==p 
            a_tr[i] = trf_yp
        else
            ip = maximum(inds)
            a_tr[i] = trf_y[ip] + ((trf_y[ip+1]-trf_y[ip]) / (trf_x[ip+1]-trf_x[ip]))*(a[i]-trf_x[ip])
        end

    end

    return a_tr
end

function invTransform(a_tr, trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
    p = length(trf_x)
    a = fill(0.0, length(a_tr))
    for i=1:length(a)
        inds = findall(a_tr[i].>trf_y)
        if length(inds)==0
            if a_tr[i] > trf_yl
                a[i] = trf_xl + ((trf_x[1]-trf_xl) / (trf_y[1]-trf_yl))*(a_tr[i]-trf_yl)
            else
                a[i] = trf_xl
            end
        elseif maximum(inds)==p
            a[i] = trf_xp
        else
            ip = maximum(inds)
            a[i] = trf_x[ip] + ((trf_x[ip+1]-trf_x[ip]) / (trf_y[ip+1]-trf_y[ip]))*(a_tr[i]-trf_y[ip])
        end
    end
    return a
end

function testTransform()
    series = rand(Uniform(0,1), 500)
    trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp = initTransform(series, 40)
    println(trf_x)
    println(trf_y)
    println(trf_xl)
    println(trf_xp)
    println(trf_yl)
    println(trf_yp)
    a = rand(Uniform(0,1), 10)
    a[1] = 0.0
    a[2]= -0.5
    a_trf = transform(a, trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
    println("a:")
    println(a)
    println("a_trf:")
    println(a_trf)
    a_trf_inv = invTransform(a_trf, trf_x, trf_y, trf_xl, trf_xp, trf_yl, trf_yp)
    println("a_trf_inv:")
    println(a_trf_inv)
    difference = a-a_trf_inv
    println(difference)
    println("Max difference: "*string(maximum(abs.(difference))))
end