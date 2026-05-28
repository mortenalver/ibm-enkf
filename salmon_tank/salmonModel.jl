using Random

mutable struct Individual
    x::Float64
    y::Float64
    v_x::Float64
    v_y::Float64
    n::Float64
    E::Float64
end

function getModelName()
    return "Salmon"
end

function setModelSettings(ms)
    # Define circular tank of diameter 8 m
    ms.xMax = 8.0
    ms.yMax = 8.0
    ms.radius = 4.0

end

function getRandomStartPos(ms)
    posOk = false
    centerPos = [ms.xMax/2.0, ms.yMax/2.0]
    pos = [ms.xMax*rand(), ms.yMax*rand()]
    while !posOk
        pos = [ms.xMax*rand(), ms.yMax*rand()]
        vFromCenter = pos - centerPos
        if sqrt(vFromCenter'*vFromCenter) < ms.radius-ms.d_wpref
            posOk = true
        end
    end
    return pos
end

function getInitialFieldLevel()
    return 0.0
end

function createIndividual(x, y, n, ms)
    ind = Individual(x, y, 0.0, 0.0, 1, 0.0)
    randomAngle = 2.0*pi*rand()
    ind.v_x = cos(randomAngle)*ms.charVel*ms.BL
    ind.v_y = sin(randomAngle)*ms.charVel*ms.BL
    return ind
end

function createIndividual(states)
    ind = Individual(states[1], states[2], states[3], states[4], states[5], states[6])
end

function step(t, dt, ind, indsArray, idx, xall, yall, X_fld, xrng, yrng, ms)

    # Update position by speed:
    ind.x = ind.x + dt*(ind.v_x)
    ind.y = ind.y + dt*(ind.v_y)

    
    # Get this inds distance from tank edge:
    cPos = [ms.xMax/2.0, ms.yMax/2.0]
    vFromC = [ind.x, ind.y] - cPos
    distFromC = sqrt(vFromC'*vFromC)
    d_w = max(0., ms.radius - distFromC)

    v_CA = [0, 0] # Unit vector towards center
    if distFromC > 0
        v_CA = -vFromC ./ sqrt(vFromC'*vFromC)
    end
    
    # Hierarchy remaining weight:
    k_remain = 1.0

    # Tank wall effect:
    v_CS = [0., 0.] # Surface
    v_CB = [0., 0.] # Bottom
    v_CW = [0., 0.] # Tank wall
    if d_w <= ms.d_wpref
        v_CW = v_CA*(ms.d_wpref - d_w)/ms.d_wpref
    end
    v_C = v_CS + v_CB + v_CW
    k_C = 1.0 # Always full weight to this component
    # Reduce remaining weigth:
    k_remain = max(0.0, k_remain - sqrt(v_C' * v_C))


    # Response towards other individuals:
    v_so_sum = [0.0, 0.0]
    if k_remain > 0.0
        xdist = abs.(xall .- ind.x)
        ydist = abs.(yall .- ind.y)
        v_so_sum = [0.0, 0.0]
        m = 0.0
        for i=1:length(xdist)
            if i==idx
                continue # Individuals do not interact with themselves
            end
            d_ij = [-xdist[i], -ydist[i]] # Vector from this individual to the other one
            d_ij_abs = sqrt(d_ij'*d_ij)  # Distance between the individuals
            if d_ij_abs <= ms.d_l # Repulsion
                v_so_sum += d_ij*(ms.d_l - d_ij_abs)
                m += 1.0
            elseif d_ij_abs <= ms.d_h # Alignment
                otherInd = indsArray[i]
                vOther = [otherInd.v_x, otherInd.v_y]
                v_so_sum += 0.5*vOther*(-d_ij_abs + ms.d_h)/(ms.d_h - ms.d_l)
                m += 1.0
            end
        end
        v_so_sum = v_so_sum / m
        #println(sqrt(v_so_sum' * v_so_sum))

    end

    # Stochastic component:
    v_ST = ms.sigma*[randn(), randn()]
    k_ST = k_remain


    # Calculate r_new:
    r_new = k_C*v_C + v_so_sum + k_ST*v_ST
    

    # Update the speed towards the newly computed speed:
    r_prev = [ind.v_x, ind.v_y]
    r_ref = ms.tau*r_prev + (1.0-ms.tau)*r_new
    spd_ref = sqrt(r_ref'*r_ref)
    r_ref_unit = r_ref
    if spd_ref > 0.0
        r_ref_unit = r_ref/sqrt(r_ref'*r_ref)
    end
    ind.v_x = ms.charVel*r_ref_unit[1]
    ind.v_y = ms.charVel*r_ref_unit[2]

    # Calculate angle between r_pref and r_ref:
    #spd_prev = sqrt(r_prev'*r_prev)
    
    #println("........")
    #rintln(r_prev)
    #println(spd_prev)
    #println(r_ref)
    #println(spd_ref)
    #println((r_prev' * r_ref)/(spd_prev*spd_ref))
    #theta = acos(max(-1.0, min(1.0, (r_prev' * r_ref)/(spd_prev*spd_ref))))
    #println(theta)
    
end


function stepAll(t, dt, indsArray, doPerturb, ms, X_fld, xrng, yrng)

    # Make x and y arrays of all individuals' positions:
    xall = zeros(Float64, length(indsArray))
    yall = zeros(Float64, length(indsArray))
    for i = 1:length(indsArray)
        xall[i] = indsArray[i].x
        yall[i] = indsArray[i].y
    end
    
    for i = 1:length(indsArray)
        ind = indsArray[i]
        step(t, dt, ind, indsArray, i, xall, yall, X_fld, xrng, yrng, ms)
    end
    return X_fld, 0
end