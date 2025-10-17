
prefix = 'C:/temp/';

dims = dlmread([prefix 'r9_resample_fieldDims.csv']);

Ubefore = dlmread([prefix 'Ubefore.csv']);
Vbefore = dlmread([prefix 'Ubefore.csv']);
Uafter = dlmread([prefix 'Uafter.csv']);
Vafter = dlmread([prefix 'Uafter.csv']);

figure, tiledlayout(2,3)
nexttile, pcolor(Ubefore'), shading flat, colorbar
nexttile, pcolor(Uafter'), shading flat, colorbar
nexttile, pcolor((Uafter-Ubefore)'), shading flat, colorbar
nexttile, pcolor(Vbefore'), shading flat, colorbar
nexttile, pcolor(Vafter'), shading flat, colorbar
nexttile, pcolor((Vafter-Vbefore)'), shading flat, colorbar


