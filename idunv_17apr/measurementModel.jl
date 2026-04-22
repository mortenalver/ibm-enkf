
function getM(dimensions, as, nPerPos)
   npos = dimensions[1]*dimensions[2]

   M = zeros(Float64, as.nmeas, npos*nPerPos)
   for i=1:as.nmeas
      state = Int(ceil(npos*rand(Float64)))
      #println(state)
      M[i, state] = 1
   end
   return M
end

function getMRegular(dims, dxy, as, nPerPos)
   npos = dims[1]*dims[2]
   stateOffsetFood = npos
   if as.speedsInStateVec
      stateOffsetFood += 2*npos
   end
   rngx = (1+trunc(Int, as.measSpacing/2)):as.measSpacing:trunc(Int, dims[1]-as.measSpacing/2)
   rngy = (1+trunc(Int, as.measSpacing/2)):as.measSpacing:trunc(Int, dims[2]-as.measSpacing/2)
   measPerPos = 1
   if as.measureFood
      measPerPos = 2
   end
   println("npos="*string(npos)*" stateOffsetFood="*string(stateOffsetFood))
   nmeas = measPerPos*length(rngx)*length(rngy)
   
   M = zeros(Float64, nmeas, npos*nPerPos)
   xmeas = zeros(Float64, nmeas)
   ymeas = zeros(Float64, nmeas)
   
   lind = LinearIndices((dims[1], dims[2]))
   idx = 0
   for i=rngx
      for j=rngy
         idx += 1
         state = lind[i,j]
         #println("i="*string(i)*", j="*string(j)*", state="*string(state))
         M[idx, state] = 1
         xmeas[idx] = i*dxy
         ymeas[idx] = j*dxy
      end
   end
   if as.measureFood
      for i=rngx
         for j=rngy
            idx += 1
            state = stateOffsetFood + lind[i,j]
            #println("i="*string(i)*", j="*string(j)*", state="*string(state))
            M[idx, state] = 1
            xmeas[idx] = i*dxy
            ymeas[idx] = j*dxy
         end
      end
   end
   return M, xmeas, ymeas
end

function getLocMatrix(dimensions, M, locDist)

   npos = dimensions[1]*dimensions[2] # Number of spatial positions in system
   sM = size(M)
   xloc = zeros(Float64, sM[2], sM[1])
   locDistSq = locDist*locDist
   # Loop over measurements:
   for i=1:sM[1]
      # Find index of this measurement:
      idx = 1
      while M[i,idx]==0
         idx = idx+1
      end
      ispatial = mod(idx, npos)
      # Find measurement coords:
      #println(idx)
      coord1 = getStateCoords(ispatial, dimensions)
      #println(coord1)
      # Loop over state variables:
      for j=1:sM[2]
         jspatial = mod(j, npos)
         coord2 = getStateCoords(jspatial, dimensions)
         distVec = [coord1[1]-coord2[1] coord1[2]-coord2[2]]
         distanceSq = distVec[1]*distVec[1] + distVec[2]*distVec[2]
         distance = sqrt(distanceSq)
         if distance < locDist
            xloc[j, i] = 1.0
         else
            xloc[j, i] = max(0.0, 1.0-(distance-locDist)/locDist)
         end

         #xloc[j, i] = exp(-distanceSq/locDistSq)

         #end
      end

   end
   return xloc
end


function getStateCoords(index, dims)
   res = zeros(Float64, 2)
   res[2] = ceil(index / dims[1])
   res[1] = index - (res[2]-1)*dims[1]
   return res
end
