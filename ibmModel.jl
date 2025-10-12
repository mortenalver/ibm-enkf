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

function getFeedConcentration(ind)
    # x* 2*y
    if ind.x>5 && ind.y>5
        return 0.5*ind.x + 1.0*ind.y
    else
        return 0.0
    end
end

function step(t, dt, ind, perturb, indsInteraction, idx, xall, yall)
    # Update position by speed with perturbations added:
    pertX = 0.0
    pertY = 0.0
    for i=1:size(perturb,1)
        # Distance from perturbation central point:
        pertVec = [perturb[i,3]-ind.x perturb[i,4]-ind.y]
        pertDist = sqrt(pertVec[1]*pertVec[1]+ pertVec[2]*pertVec[2])
        distFactor = 1.0/(1.0 + exp(pertDist - 1.0))
        pertX = pertX + perturb[i,1]*distFactor
        pertY = pertY + perturb[i,2]*distFactor
    end
    ind.x = ind.x + dt*(ind.v_x + pertX)
    ind.y = ind.y + dt*(ind.v_y + pertY)

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
    # Standard speed:
    v_norm = ind.normSpeed
    # Find speed vector towards setpoint:
    v_x_new = 0
    v_y_new = 0
    if distToRef > 0
        v_x_new = v_norm*vecToRef[1]/distToRef
        v_y_new = v_norm*vecToRef[2]/distToRef
    end
    # Update the speed towards the newly computed speed:
    ind.v_x += dt*(v_x_new - ind.v_x)
    ind.v_y += dt*(v_y_new - ind.v_y)

    # If individual interaction is activated, orient away from close individuals:
    if indsInteraction
        distThresh = 0.2
        dxInt = 0.0
        dyInt = 0.0
        xdist = abs.(xall .- ind.x)
        ydist = abs.(yall .- ind.y)
        for i=1:length(xdist)
            if i==idx
                continue # Individuals do not interact with themselves
            end
            # Check if x+y distance is less than a threshold: 
            if (xdist[i]+ydist[i]) < (1.5*distThresh)
                dist = sqrt(xdist[i]*xdist[i] + ydist[i]*ydist[i])
                if dist < distThresh
                    # Add a little to x and y speeds to move away from the close individual:
                    abxfac = abs((0.5*distThresh)/xdist[i])
                    abyfac = abs((0.5*distThresh)/ydist[i])
                    
                    dxInt = dxInt + 0.1*sign(ind.x - xall[i])*min(1.0, abxfac)
                    dyInt = dyInt + 0.1*sign(ind.y - yall[i])*min(1.0, abyfac)
                end
            end
        end

        # Add the summed adjustments to speed:
        ind.v_x += dt*dxInt
        ind.v_y += dt*dyInt
        
    end

    # Update energy level:
    X = getFeedConcentration(ind)
    f = X/(X + 10)
    ind.E = ind.E + dt*(f - 0.25*ind.E)

end

function stepAll(t, dt, indsArray, perturb, ms)

    xall = []
    yall = []
    if ms.indsInteraction
        # Make x and y arrays of all individuals' positions:
        xall = zeros(Float64, length(indsArray))
        yall = zeros(Float64, length(indsArray))
        for i = 1:length(indsArray)
            xall[i] = indsArray[i].x
            yall[i] = indsArray[i].y
        end
    end

    for i = 1:length(indsArray)
        ind = indsArray[i]
        #xind = round(Int64, ind.x)
        #yind = round(Int64, ind.y)
        #println(i)
        #if i==1
        #    println(ufield[xind, yind])
        #end
        step(t, dt, ind, perturb, ms.indsInteraction, i, xall, yall)
        
    end
    
end