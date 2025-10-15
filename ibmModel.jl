using Random

mutable struct Individual
    x::Float64
    y::Float64
    v_x::Float64
    v_y::Float64
    E::Float64
    normSpeed::Float64
    n::Float64
    refX::Float64
    refY::Float64
    lastArea::Int
end

function createIndividual(x, y, normSpeed, n)
    ind = Individual(x, y, 0.0, 0.0, 0.0, normSpeed, n, 0.0, 0.0, -1)
    return ind
end


function getSpeedMult(X)
    # Get speed multiplier as a function of local food concentration:
    # X is in the interval [0, 1], and speed decreases for high values.
    return 1.3 - 0.6*X


end

function step(t, dt, ind, perturb, idx, xall, yall, X_fld, xrng, yrng, ms)
    
    # Get local feed concentration:
    # First find index in X field:
    dxy = xrng[2]-xrng[1]
    ix = max(1, min(Int(floor((xall[idx] - xrng[1])/dxy)), size(X_fld,1)))
    iy = max(1, min(Int(floor((yall[idx] - yrng[1])/dxy)), size(X_fld,2)))
    X = X_fld[ix,iy]
    
    # Get food-dependent speed multiplier:
    XSpeedMult = getSpeedMult(X)

    # Update position by speed:
    ind.x = ind.x + dt*(ind.v_x)
    ind.y = ind.y + dt*(ind.v_y)

    # Updated speed component setpoints, to be computed:
    v_x_new = 0.0
    v_y_new = 0.0

    if ms.migration
        # Figure out where the setpoint should be:
        intervalLength = 7 # Time interval in which we aim towards one of the target areas before switching
        targetArea = 1+mod(floor(t/intervalLength),4) # Find target area by the current time
        # Check if the target area has changed. If so, find new set point:
        if targetArea != ind.lastArea
            ind.lastArea = targetArea
            # Determine the target rectangle:
            targetRect = [13 18; 8 13] # Area 1 (x1 x2; y1 y2)
            if targetArea==2
                targetRect = [13 18; 2 7] # Area 2
            elseif targetArea==3
                targetRect = [2 7; 2 7] # Area 2
            elseif targetArea==4
                targetRect = [2 7; 8 13] # Area 2
            end        
            # Choose a random setpoint in the target 
            #ind.refX = targetRect[1,1] + (targetRect[1,2]-targetRect[1,1])*rand(Float64)
            #ind.refY = targetRect[2,1] + (targetRect[2,2]-targetRect[2,1])*rand(Float64)
            ind.refX = 0.5*(targetRect[1,2]+targetRect[1,1]) + 0.25*(targetRect[1,2]-targetRect[1,1])*randn(Float64)
            ind.refY = 0.5*(targetRect[2,2]+targetRect[2,1]) + 0.25*(targetRect[2,2]-targetRect[2,1])*randn(Float64)
        end
        # Find vector from current position to setpoint:
        vecToRef = [ind.refX-ind.x; ind.refY-ind.y]
        # Find distance to setpoint:
        distToRef = sqrt(vecToRef[1]*vecToRef[1] + vecToRef[2]*vecToRef[2])
        # Adjusted standard speed:
        v_norm = XSpeedMult*ind.normSpeed
        # Find speed vector towards setpoint:
        v_x_new = 0
        v_y_new = 0
        if distToRef > 0
            v_x_new = v_norm*vecToRef[1]/distToRef
            v_y_new = v_norm*vecToRef[2]/distToRef
        end
    else # No migration
        v_x_new = randn()
        v_y_new = randn()
    end

    # Add random perturbations to the speed components:
    pertX = 0.0
    pertY = 0.0
    for i=1:size(perturb,1)
        # Distance from perturbation central point:
        pertVec = [perturb[i,3]-ind.x perturb[i,4]-ind.y]
        pertDist = sqrt(pertVec[1]*pertVec[1]+ pertVec[2]*pertVec[2])
        pertDist = pertDist / 5.0
        distFactor = exp(-(pertDist*pertDist))
        pertX = pertX + perturb[i,1]*distFactor
        pertY = pertY + perturb[i,2]*distFactor
    end
    v_x_new += pertX
    v_y_new += pertY

    # Update the speed towards the newly computed speed:
    ind.v_x += dt*0.6*(v_x_new - ind.v_x)
    ind.v_y += dt*0.6*(v_y_new - ind.v_y)

    # If individual interaction is activated, orient away from close individuals:
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
            if (xdist[i]+ydist[i]) < (1.5*ms.indsInteractionThresh)
                dist = sqrt(xdist[i]*xdist[i] + ydist[i]*ydist[i])
                if dist < ms.indsInteractionThresh
                    # Add a little to x and y speeds to move away from the close individual:
                    abxfac = abs((0.5*ms.indsInteractionThresh)/xdist[i])
                    abyfac = abs((0.5*ms.indsInteractionThresh)/ydist[i])
                    
                    dxInt = dxInt + ms.indsInteractionStrength*sign(ind.x - xall[i])*min(1.0, abxfac)
                    dyInt = dyInt + ms.indsInteractionStrength*sign(ind.y - yall[i])*min(1.0, abyfac)
                end
            end
        end

        # Add the summed adjustments to speed:
        ind.v_x += dt*dxInt
        ind.v_y += dt*dyInt
        
    end

    # Update energy level:
    f = X/(X + 0.5)
    ind.E = ind.E + dt*(f - 0.25*ind.E)
    X_fld[ix, iy] = max(0.0, X_fld[ix, iy]-dt*0.01*f)

end

function stepAll(t, dt, indsArray, perturb, ms, X_fld, xrng, yrng)

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

    for i = 1:length(indsArray)
        ind = indsArray[i]
        #xind = round(Int64, ind.x)
        #yind = round(Int64, ind.y)
        #println(i)
        #if i==1
        #    println(ufield[xind, yind])
        #end
        step(t, dt, ind, perturb, i, xall, yall, X_fld, xrng, yrng, ms)
        
    end
   
    # Perturb feed concentration field:
    X_pert = getRandomField([size(X_fld,1) size(X_fld,2)], 0.2, 25, 6)
    X_fld = X_fld + dt.*X_pert

    # Add food everywhere to a maximum of 2:
    X_fld = X_fld + fill(dt*0.025, size(X_fld,1), size(X_fld,2))
    X_fld = min.(X_fld, 2.0)
    return X_fld
end