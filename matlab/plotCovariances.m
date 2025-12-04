
% Read ensemble:
%X = readmatrix('C:/temp/X_f_all.csv');
X = readmatrix('D:/work/ibm-enkf/test1/X_f_all.csv');
N = size(X,2);

% Read dimensions:
dims_and_int = readmatrix("C:/temp/d_indi2500_12_resample_fieldDims.csv");
dims = dims_and_int(1:2);
npos = dims(1)*dims(2);

figure, tiledlayout('flow','TileSpacing','compact');

for ip = 5:25:npos
    %point = [10 10];
    %ip = sub2ind(dims, point(1), point(2));
    [ii,jj] = ind2sub(dims, ip);
    s1 = X(ip,:);

    Cx = zeros(dims);
    
    for i=1:dims(1)
        for j=1:dims(2)
            jp = sub2ind(dims, i, j);
            s2 = X(jp,:);
            cc = corrcoef([s1' s2']);
            if ~isnan(cc(1,2))
                Cx(i,j) = cc(1,2);
            end
        end
    end
    
    nexttile,pcolor(Cx'), shading flat
    title(num2str(ii)+" / "+num2str(jj));

end