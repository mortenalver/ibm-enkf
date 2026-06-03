using Random

mutable struct Individual
    x::Float64
    y::Float64
    v_x::Float64
    v_y::Float64
    E::Float64
    n::Float64
end

function getIndStateVec(ind)
    return [ind.x, ind.y, ind.v_x, ind.v_y, ind.n, ind.E]
end


function createIndividual(x, y, n, ms)
    ind = Individual(x, y, 0.0, 0.0, 0.0, n)
    return ind
end

function createIndividual(vec)
    ind = Individual(vec[1], vec[2], vec[3], vec[4], vec[5], vec[6])
    return ind
end

function setModelSettings(ms)
    ms.xMax = 20
    ms.yMax = 15
    ms.k_X = 2.5 # Strength of motion towards food gradient.
    ms.X_a = 0.01 # Global addition rate of food.
    ms.indsInteraction = true # If true, individuals will be repulsed from each other at close distances (not very optimized, so makes model slower)
    ms.indsAlignThresh = 0.3 # Distance threshold for individual alignment.
    ms.indsRepulseThresh = 0.08 # Distance threshold for individual repulsion.
    ms.indsAttractThresh = 0.75 # Distance threshold for individual attraction.
    ms.n_wuv = 4 # Number of perturbations for u/v
    ms.n_wX = 2 # Number of perturbations for X
    ms.sigma_uv = 1.0 # 3.0 # Intensity of u/v random perturbations
    ms.sigma_X = 0.8 # Intensity of X random perturbations
    ms.d_wuv = 5.0
    ms.d_wX = 6.0
    ms.indsAlignStrength = 0.01 # Strength of individual interaction
    ms.indsRepulseStrength = 0.4#0.08 # Strength of individual interaction
    ms.indsAttractStrength = 0.002 # Strength of individual interaction
    ms.pullTowardsCOG = false # If true, individuals will be pulled towards the center of gravity of the population
    ms.pullTowardsCOGStrength = 0.75#0.3 # Strength of the pull towards COG if activated
    ms.speedUpdateRate = 0.4 # 0.6 # Multiplier for speed update - lower means more intertia in speed updates
    ms.nInd = 5000 #6000 # Number of individuals
    ms.nPerInd = 1.0 # individuals per super individual

end

function initializeSettings(ms, dims)

end

function getModelName()
    return "Test"
end

function getRandomStartPos(ms)
    return [3.1 + 2*randn(Float64), 9.1 + 2*randn(Float64)]
end

function getInitialFieldLevel()
    return 1.0
end


function createIndividual(states)
    ind = Individual(states[1], states[2], states[3], states[4], states[5], states[6])
end

function getSpeedMult(X)
    # Get speed multiplier as a function of local food concentration:
    # X is in the interval [0, 1], and speed decreases for high values.
    return 1.3 - 0.6*X


end

function step(t, dt, ind, indsArray, perturb, idx, xall, yall, X_fld, xrng, yrng, cog, ms)
    
    # Get local feed concentration:
    # First find index in X field:
    dxy = xrng[2]-xrng[1]
    ix = max(1, min(Int(floor((xall[idx] - xrng[1])/dxy)), size(X_fld,1)))
    iy = max(1, min(Int(floor((yall[idx] - yrng[1])/dxy)), size(X_fld,2)))
    X = X_fld[ix,iy]
    
    # Update position by speed:
    ind.x = ind.x + dt*(ind.v_x)
    ind.y = ind.y + dt*(ind.v_y)

    # Updated speed component setpoints, to be computed:
    v_x_new = 0.0
    v_y_new = 0.0

    # We first decide on a vector tendency, then add random movement to it.
    tend_x = 0.0
    tend_y = 0.0
    # If we are near an edge, add a tendency away from the edge:
    bnd = 3
    k_bnd = 1.0
    if ind.x < bnd
        tend_x += k_bnd*(bnd-ind.x)/bnd
    end
    if ms.xMax-ind.x < bnd
        tend_x -= k_bnd*(bnd-(ms.xMax-ind.x))/bnd
    end
    if ind.y < bnd
        tend_y += k_bnd*(bnd-ind.y)/bnd
    end
    if ms.yMax-ind.y < bnd
        tend_y -= k_bnd*(bnd-(ms.yMax-ind.y))/bnd
    end
    # Then, add a tendency depending along the local X_fld gradient:
    if ix>1 && ix<size(X_fld,1)
        tend_x += ms.k_X*(X_fld[ix+1,iy]-X_fld[ix-1,iy])
    end
    if iy>1 && iy<size(X_fld,2)
        tend_y += ms.k_X*(X_fld[ix,iy+1]-X_fld[ix,iy-1])
    end

    # Then, if activated, add a tendency towards the center of the population:
    if ms.pullTowardsCOG
        vecToCOG = [cog[1]-ind.x cog[2]-ind.y]
        lVec = sqrt(vecToCOG[1]*vecToCOG[1] + vecToCOG[2]*vecToCOG[2])
        # If we are not exactly at the COG, add a pull in that direction:
        if maximum(lVec) > 9
            normVec = vecToCOG/lVec
            tend_x += ms.pullTowardsCOGStrength*normVec[1]
            tend_y += ms.pullTowardsCOGStrength*normVec[2]
        end
    end

    v_x_new = tend_x + 0.5*randn()
    v_y_new = tend_y + 0.5*randn()


    # Add random perturbations to the speed components:
    pertX = 0.0
    pertY = 0.0
    for i=1:size(perturb,1)
        # Distance from perturbation central point:
        pertVec = [perturb[i,3]-ind.x perturb[i,4]-ind.y]
        pertDist = sqrt(pertVec[1]*pertVec[1]+ pertVec[2]*pertVec[2])
        pertDist = pertDist / ms.d_wuv #/ 10.0 #5.0
        distFactor = exp(-(pertDist*pertDist))
        pertX = pertX + perturb[i,1]*distFactor
        pertY = pertY + perturb[i,2]*distFactor
    end
    v_x_new += pertX
    v_y_new += pertY

    
    # If individual interaction is activated, align with close individuals:
    totInteractions = 0
        
    if ms.indsInteraction
        dxInt = 0.0
        dyInt = 0.0
        xdist = abs.(xall .- ind.x)
        ydist = abs.(yall .- ind.y)
        for i=1:length(xdist)
            if i==idx
                continue # Individuals do not interact with themselves
            end
            # Check if x+y distance is less than a threshold: 
            if (xdist[i]+ydist[i]) < (ms.indsAttractThresh)
                dist = sqrt(xdist[i]*xdist[i] + ydist[i]*ydist[i])
                if dist < ms.indsAlignThresh
                    totInteractions = totInteractions + 1

                    # If closer than the repulse threshold, add a little to x and y speeds 
                    # to move away from the close individual:
                    if dist < ms.indsRepulseThresh
                        #abxfac = abs(ms.indsRepulseThresh/xdist[i])
                        #abyfac = abs(ms.indsRepulseThresh/ydist[i])
                        abfac = (ms.indsRepulseThresh - dist)/ms.indsRepulseThresh
                        dxInt = dxInt + ms.indsRepulseStrength*sign(ind.x - xall[i])*abfac*xdist[i]/dist
                        dyInt = dyInt + ms.indsRepulseStrength*sign(ind.y - yall[i])*abfac*ydist[i]/dist


                    end
                    
                    # Add a little to x and y speeds to align speed with the close individual:
                    relDist = dist/ms.indsAlignThresh
                    nearInd = indsArray[i]
                    dxInt += ms.indsAlignStrength*(1.0-relDist)*(nearInd.v_x - ind.v_x)
                    dyInt += ms.indsAlignStrength*(1.0-relDist)*(nearInd.v_y - ind.v_y)
                    

                elseif dist < ms.indsAttractThresh
                    totInteractions = totInteractions + 1

                    # If further from align threshold but closer than attract threshold, add a little to
                    # x and y speeds to move towards the close individual:
                    dxInt = dxInt + ms.indsAttractStrength*sign(xall[i] - ind.x)*xdist[i]/dist
                    dyInt = dyInt + ms.indsAttractStrength*sign(yall[i] - ind.y)*ydist[i]/dist
                    #dxInt = dxInt + ms.indsAttractStrength*((ms.indsAttractThresh-dist)/(ms.indsAttractThresh-ms.indsAlignThresh))*sign(xall[i] - ind.x)*xdist[i]/dist
                    #dyInt = dyInt + ms.indsAttractStrength*((ms.indsAttractThresh-dist)/(ms.indsAttractThresh-ms.indsAlignThresh))*sign(yall[i] - ind.y)*ydist[i]/dist
                end
            end
        end
        
        # Add the summed adjustments to the new speed speed:
        v_x_new += dxInt
        v_y_new += dyInt

    end

    
    # Update the speed towards the newly computed speed:
    ind.v_x += dt*ms.speedUpdateRate*(v_x_new - ind.v_x)
    ind.v_y += dt*ms.speedUpdateRate*(v_y_new - ind.v_y)
    

    # Update energy level:
    f = X/(X + 0.5)
    ind.E = ind.E + dt*(f - 0.25*ind.E)
    #println(string(ind.E)*" / "*string(f))
        
    X_fld[ix, iy] = max(0.0, X_fld[ix, iy]-dt*0.015*f)

    return totInteractions
end

function stepAll(t, dt, indsArray, doPerturb, ms, X_fld, xrng, yrng)

    xlim = [0, ms.xMax]
    ylim = [0, ms.yMax]
    
    # Roll random numbers used to perturb individuals' speeds:
    perturb = zeros(Float64, ms.n_wuv, 4)
    if doPerturb
        for ptI = 1:size(perturb,1)
            perturb[ptI,:] = [ms.sigma_uv*randn(Float64), ms.sigma_uv*randn(Float64),
                xlim[1]+rand(Float64)*(xlim[2]-xlim[1]), ylim[1]+rand(Float64)*(ylim[2]-ylim[1])]
        end
    end

    xall = []
    yall = []
    #if ms.indsInteraction
    # Make x and y arrays of all individuals' positions:
    xall = zeros(Float64, length(indsArray))
    yall = zeros(Float64, length(indsArray))
    for i = 1:length(indsArray)
        xall[i] = indsArray[i].x
        yall[i] = indsArray[i].y
    end
    #end

    # If necessary, compute the center of gravity:
    cog = [0.0 0.0]
    if ms.pullTowardsCOG
        xsum = 0.0
        ysum = 0.0
        ncount = 0.0
        for i in eachindex(indsArray)
            ind = indsArray[i]
            xsum += ind.n*ind.x
            ysum += ind.n*ind.y
            ncount += ind.n
        end
        cog[1] = xsum/ncount
        cog[2] = ysum/ncount
    end

    totInteractions = 0
    for i = 1:length(indsArray)
        ind = indsArray[i]

        interactions = step(t, dt, ind, indsArray, perturb, i, xall, yall, X_fld, xrng, yrng, cog, ms)
        totInteractions += interactions
    end
    #println("Average interactions per ind: "*string(totInteractions/length(indsArray)))

    # Perturb feed concentration field:
    X_pert = getRandomField([size(X_fld,1) size(X_fld,2)], ms.sigma_X, ms.n_wX, ms.d_wX, 0.0)
    X_fld = X_fld + dt.*X_pert 

    # Add food everywhere to a maximum of 4:
    X_fld = X_fld + fill(dt*ms.X_a, size(X_fld,1), size(X_fld,2))
    X_fld = min.(X_fld, 4.0)
    return X_fld, totInteractions/length(indsArray)
end
