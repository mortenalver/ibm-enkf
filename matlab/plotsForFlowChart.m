prefix = 'C:/temp/';

X_pre = readmatrix(prefix+"dens_before.csv");
X_post = readmatrix(prefix+"dens_after.csv");

x_e1 = dlmread(prefix+"forplotting_resample_e1X.csv");
y_e1 = dlmread(prefix+"forplotting_resample_e1Y.csv");

x_twin = dlmread(prefix+"forplotting_resample_twinX.csv");
y_twin = dlmread(prefix+"forplotting_resample_twinY.csv");
dens_twin = dlmread(prefix+"forplotting_resample_twinDens.csv");

figure,tiledlayout('flow')
dx = 0.5;
idx = 40; % Time step to plot
lims = [0 50];
xlims = [dx 12];
ylims = [dx 15+dx]


nexttile, scatter(x_e1(:,idx-1), y_e1(:,idx-1), 0.5), axis image, xlim([0 20]), ylim([0 15]), title('E1 pre')
xlim(xlims), ylim(ylims)
nexttile, scatter(x_e1(:,idx), y_e1(:,idx), 0.5), axis image, xlim([0 20]), ylim([0 15]), title('E1 post')
xlim(xlims), ylim(ylims)
nexttile, pcolor(dx*(1:size(X_pre,1)), dx*(1:size(X_pre,2)), X_pre'), axis image, shading flat, clim(lims)
xlim(xlims), ylim(ylims)
nexttile, pcolor(dx*(1:size(X_pre,1)), dx*(1:size(X_pre,2)), X_post'), axis image, shading flat, clim(lims)
xlim(xlims), ylim(ylims)
nexttile, scatter(x_twin(:,idx), y_twin(:,idx), 0.5), axis image, xlim([0 20]), ylim([0 15]), title('Twin')
xlim(xlims), ylim(ylims)
nexttile
dFieldT = reshape(dens_twin(:,idx), dims(1), dims(2));
pcolor(dx*(1:size(X_pre,1)), dx*(1:size(X_pre,2)), dFieldT'), shading flat, axis image
clim(lims)
title('Density field (twin)')
xlim(xlims), ylim(ylims)