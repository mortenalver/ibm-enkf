function X = downScale(M, nsub)
% Downscale matrix M by a factor nsub by averaging nsub*nsub squares.

X = zeros(ceil(size(M,1)/nsub), ceil(size(M,2)/nsub));

for i=1:size(X,1)
    for j=1:size(X,2)
        square = M(1+nsub*(i-1):min(size(M,1), nsub*i), ...
            1+nsub*(j-1):min(size(M,2), nsub*j));
        X(i,j) = mean(square(:));
    end
end
