function rmsVal = myRms(M, dims, dsFactor)

if nargin < 3
    dsFactor = 1;
end

matrix = reshape(M, dims(1), dims(2));
mdown = downScale(matrix, dsFactor);
rmsVal = rms(mdown(:));