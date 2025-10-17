
ot = dlmread('C:/temp/OT.csv');
orig = dlmread('C:/temp/orig.csv');
orig0 = dlmread('C:/temp/orig0.csv');
corre = dlmread('C:/temp/corr.csv');
corre0 = dlmread('C:/temp/corr0.csv');
postdens = dlmread('C:/temp/postdens.csv');

figure
nexttile, pcolor(orig'), shading flat, colorbar, title('Original subsetted dens'), clim([0 10])
nexttile, pcolor(corre'), shading flat, colorbar, title('Corrected subsetted dens'), clim([0 10])
nexttile, pcolor(orig0'), shading flat, colorbar, title('Original dens'), clim([0 1])
nexttile, pcolor(corre0'), shading flat, colorbar, title('Corrected dens'), clim([0 1])
nexttile, pcolor((corre-orig)'), shading flat, colorbar, title('Corrected - Orig'), clim([-0.1 0.1])

fromCoord = [2 2]
indCoord = sub2ind(size(orig), fromCoord(1), fromCoord(2))
orig_corre_point = [orig(fromCoord(1), fromCoord(2)) corre(fromCoord(1), fromCoord(2))]

toVal = ot(indCoord,:)

nexttile, plot(toVal)
%nexttile, plot(ot(:,indCoord))

[m, I] = maxk(toVal,2)
[secondMosttoCoordI, secondMosttoCoordJ] = ind2sub(size(orig), I(2))

% for i=1:size(ot,2)
% 
% 
% end


%%
% Go through OT matrix and distribute into new matrix:
testfield = zeros(size(corre));
for i=1:size(ot,1)
    for j=1:size(ot,2)
        testfield(j) = testfield(j) + ot(i,j);
    end
end

nexttile,pcolor(testfield'),shading flat, colorbar, title('Recreated from OT')
nexttile,pcolor((corre-testfield)'),shading flat, colorbar, title('Deviation in recreation')
nexttile,pcolor(postdens'),shading flat, colorbar, title('Dens after OT')
