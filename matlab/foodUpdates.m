
prefix = 'C:/temp/';

dims = dlmread([prefix 'r9_resample_fieldDims.csv']);

xloc = dlmread([prefix 'Xloc.csv']);

Xbefore = dlmread([prefix 'Xbefore.csv']);
Xafter = dlmread([prefix 'Xafter.csv']);

figure, tiledlayout(2,2)
nexttile, pcolor(Xbefore'), shading flat, colorbar, axis image
nexttile, pcolor(Xafter'), shading flat, colorbar, axis image
nexttile, pcolor((Xafter-Xbefore)'), shading flat, colorbar, axis image


