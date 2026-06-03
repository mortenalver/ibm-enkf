using Random

mutable struct Individual
    x::Float64
    y::Float64
    v_x::Float64
    v_y::Float64
    n::Float64
    E::Float64
    mode::Int
end

function getModelName()
    return "Salmon"
end

function getIndStateVec(ind)
    return [ind.x, ind.y, ind.v_x, ind.v_y, ind.n, ind.E, ind.mode]
end


function createIndividual(x, y, n, ms)
    ind = Individual(x, y, 0.0, 0.0, 1, 0.0, 0)
    randomAngle = 2.0*pi*rand()
    ind.v_x = cos(randomAngle)*ms.charVel*ms.BL
    ind.v_y = sin(randomAngle)*ms.charVel*ms.BL
    return ind
end

function createIndividual(vec)
    mode = Int(round(vec[7]))
    ind = Individual(vec[1], vec[2], vec[3], vec[4], vec[5], vec[6], mode)
    return ind
end

function setModelSettings(ms)

end

function setStateVector(ind, vec)
    ind.x = vec[1]
    ind.y = vec[2]
    ind.v_x = vec[3]
    ind.v_y = vec[4]
    ind.n = vec[5]
    ind.E = vec[6]
    ind.mode = Int(round(vec[7]))
    #println(string(ind.x)*"  "*string(ind.mode))
end

function initializeSettings(ms, dims)
    global feedingField = fill(0.0, dims[1], dims[2])
    feedRad = 1.1
    cPos = [ms.xMax/2.0, ms.yMax/2.0]
    nFeedPoints = 0
    for i=1:dims[1]
        for j=1:dims[2]
            vFromC = ms.dxy.*([i, j].-0.5) - cPos
            distFromC = sqrt(vFromC' * vFromC)
            if distFromC <= feedRad
                feedingField[i,j] = 1.0
                nFeedPoints += 1
            end
        end
    end
    feedingField /= nFeedPoints
    
end

function getRandomStartPos(ms)
    posOk = false
    centerPos = [ms.xMax/2.0, ms.yMax/2.0]
    pos = [ms.xMax*rand(), ms.yMax*rand()]
    while !posOk
        pos = [ms.xMax*rand(), ms.yMax*rand()]
        vFromCenter = pos - centerPos
        if sqrt(vFromCenter'*vFromCenter) < ms.radius-1.0*ms.d_wpref
            posOk = true
        end
    end
    return pos
end

function getInitialFieldLevel()
    return 0.0
end


# function createIndividual(states)
#     ind = Individual(states[1], states[2], states[3], states[4], states[5], states[6])
# end

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
    k_C = 2.0 # Always full weight to this component
    if d_w < 0.25*ms.d_wpref
        k_C = 5.0 # Extra weight when really close to the wall.
    end
    # Reduce remaining weigth:
    k_remain = max(0.0, k_remain - sqrt(v_C' * v_C))


    # Feeding behaviour state machine:
    maxFood = maximum(X_fld)
    feedPresent = maxFood > 1e-3 # True if there is feed in the system
    if ind.mode == 0     # normal
        # Check if there is feed in the system:
        if feedPresent
            # FOR NOW, ASSUMING PROBABILITY INDEPENDENT OF DISTANCE TO FEEDING
            # LOCATION SINCE WE ARE IN A TANK RATHER THAN A CAGE
            transitionProb = dt*0.1
            if rand() < transitionProb
                ind.mode = 1 # Transition to approach
            end

        end
    elseif ind.mode == 1 # approach
        if !feedPresent
            # No feed, så gå back to normal model
            ind.model = 0
        end
    elseif ind.mode == 2 # manipulate
        
    elseif ind.mode == 3 # satiated

    end


    # Response towards food:
    v_F = [0, 0]
    if ind.mode == 0     # normal
        # No food effect
    elseif ind.mode == 1 # approach
        # Go towards food maximum:
        indMaxFood = argmax(X_fld)
        posMaxFood = ms.dxy.*([indMaxFood[1], indMaxFood[2]].-0.5)
        v_F = posMaxFood-[ind.x, ind.y]
        if sqrt(v_F' * v_F) > 0
            v_F /= sqrt(v_F' * v_F)
        end

        # Check if the fish succeeds in eating.
        # First find the local food concentration:
        ix = max(1, min(Int(floor((xall[idx] - xrng[1])/ms.dxy)), size(X_fld,1)))
        iy = max(1, min(Int(floor((yall[idx] - yrng[1])/ms.dxy)), size(X_fld,2)))
        X = X_fld[ix,iy]
        if X > ms.pelletWeight
            captureProb = dt*0.5
            if rand() < captureProb
                # Captured pellet
                ind.E += ms.pelletWeight
                X_fld[ix,iy] -= ms.pelletWeight
                ind.mode = 2 # State transition to manipulate
            end
            
        end
    elseif ind.mode == 2 # manipulate

    elseif ind.mode == 3 # satiated

    end

    # Response towards other individuals:
    d_h = ms.d_h*ms.BL
    d_l = ms.d_l*ms.BL
    #println(d_h)
    #println(d_l)
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
            if d_ij_abs <= d_l # Repulsion
                v_so_sum += d_ij*(d_l - d_ij_abs)
                m += 1.0
                if any(isnan.(v_so_sum))
                    println("1: ", v_so_sum)
                end
            elseif d_ij_abs <= d_h # Alignment
                otherInd = indsArray[i]
                vOther = [otherInd.v_x, otherInd.v_y]
                v_so_sum += 0.5*vOther*(-d_ij_abs + d_h)/(d_h - d_l)
                m += 1.0
                if any(isnan.(v_so_sum))
                    println("2: ", v_so_sum)
                end
            end
        end
        if m>0
            v_so_sum = v_so_sum / m
        end
        l_v_so = sqrt(v_so_sum' * v_so_sum)
        if any(isnan.(v_so_sum))
            println("3: ", v_so_sum)
        end
        if l_v_so > k_remain
            v_so_sum *= k_remain/l_v_so
            k_remain = 0.0
        else
            k_remain -= l_v_so
        end
    end

    

    # Stochastic component:
    v_ST = ms.sigma*[randn(), randn()]
    k_ST = k_remain


    # Calculate r_new:
    r_new = k_C*v_C + v_F + v_so_sum + k_ST*v_ST
    

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
    
    # Update all individuals:
    for i = 1:length(indsArray)
        ind = indsArray[i]
        step(t, dt, ind, indsArray, i, xall, yall, X_fld, xrng, yrng, ms)
    end

    # Add food:
    if t > ms.feedingStart && ms.feedingRate > 0.0
        for i1 = 1:size(feedingField,1)
            for i2 = 1:size(feedingField,2)
                X_fld[i1,i2] += dt*feedingField[i1,i2]
            end
        end

    end

    return X_fld, 0
end