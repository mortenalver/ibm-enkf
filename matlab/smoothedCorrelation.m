function corre = smoothedCorrelation(f1, f2, dims, doplot)

fld1 = reshape(f1, dims);
fld2 = reshape(f2, dims);

if doplot
    figure
    nexttile, mesh(fld1)
    nexttile, mesh(fld2)
end
winsize = 3;
fld1 = smoothdata2(fld1,"gaussian",winsize);
fld2 = smoothdata2(fld2,"gaussian",winsize);
if doplot
    nexttile, mesh(fld1)
    nexttile, mesh(fld2)
end

corrval = corrcoef([reshape(fld1,numel(fld1),1) reshape(fld2,numel(fld2),1)]);

if ~isnan(corrval(1,2))
    corre = corrval(1,2);
else
    corre = 0;
end

% if doplot
%     1+1
% end
