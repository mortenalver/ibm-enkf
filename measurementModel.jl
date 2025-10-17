
function getM(dimensions, as)
   nmeas = 2*800
   npos = dimensions[1]*dimensions[2]
   nPerPos = 1 
   if as.speedsInStateVec
      nPerPos = 3
   end
   M = zeros(Float64, nmeas, npos*nPerPos)
   for i=1:nmeas
      state = Int(ceil(npos*rand(Float64)))
      #println(state)
      M[i, state] = 1
   end
   return M
end

function getLocMatrix(dimensions, M, locDist)

   sM = size(M)
   xloc = zeros(Float64, sM[2], sM[1])

   # Loop over measurements:
   for i=1:sM[1]
      # Find index of this measurement:
      idx = 1
      while M[i,idx]==0
         idx = idx+1
      end
      # Find measurement coords:
      #println(idx)
      coord1 = getStateCoords(idx, dimensions)
      #println(coord1)
      # Loop over state variables:
      for j=1:sM[2]
         coord2 = getStateCoords(j, dimensions)
         distVec = [coord1[1]-coord2[1] coord1[2]-coord2[2]]
         distance = sqrt(distVec[1]*distVec[1] + distVec[2]*distVec[2])
         if distance < locDist
            xloc[j, i] = 1.0
         end
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
