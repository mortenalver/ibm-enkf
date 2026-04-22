
%prefix = 'D:/work/ibm-enkf/test1/perturb_5_';
prefix = 'D:/work/ibm-enkf/test1/run_2026_3_resample_';
%prefix = 'C:/temp/perturb_1_resample_';
%prefix = 'C:/temp/d_gtwin_perturb_1_resample_';
%prefix = 'C:/temp/d_gtwin_indi2500_13_resample_';
%prefix = 'C:/temp/indi2500_10_resample_';
%prefix = 'C:/temp/nomigr11_resample_';
%prefix = 'C:/temp/d_nomigr_long1_resample_';

animate = 1;
plotInd = 25;

plotDensField = 1;
cols = 1;
figW = 650;
if plotDensField
    cols = 2;
    figW = 1200;
end

dims = dlmread([prefix 'fieldDims.csv']);

% Read twin values:

dt = 1*0.2;
if length(dims) > 3
    dt = dims(4);
end
speedup = 1;

x_twin = dlmread([prefix 'twinX.csv']);
y_twin = dlmread([prefix 'twinY.csv']);
E_twin = dlmread([prefix 'twinE.csv']);
N_twin = dlmread([prefix 'twinN.csv']);
U_twin = dlmread([prefix 'twinU.csv']);
V_twin = dlmread([prefix 'twinV.csv']);
food_twin = dlmread([prefix 'twinFood.csv']);
dens_twin = dlmread([prefix 'twinDens.csv']);
energy_twin = dlmread([prefix 'twinEnergy.csv']);
Xfld_twin = dlmread([prefix 'twinXfld.csv']);

enkfField = dlmread([prefix 'enkfField.csv']);
% Read ensemble values:
x_1 = dlmread([prefix 'e1X.csv']);
y_1 = dlmread([prefix 'e1Y.csv']);
E_1 = dlmread([prefix 'e1E.csv']);
N_1 = dlmread([prefix 'e1N.csv']);
food_1 = dlmread([prefix 'e1Food.csv']);
dens_e = dlmread([prefix 'eDens.csv']);
densStd_e = dlmread([prefix 'eDensStd.csv']);
energy_e = dlmread([prefix 'eEnergy.csv']);
Xfld_e = dlmread([prefix 'eXfld.csv']);


%%


v = VideoWriter('anim','MPEG-4');
v.FrameRate = 4; % Default 30
open(v);

figure('Renderer', 'painters', 'Position', [10 50 figW 450])
ncol = 4;
tiledlayout(1,cols)

if animate > 0
    range = 1:speedup:size(x_twin,2)
else
    range = plotInd;
end

rmsValues = zeros(length(range),4);
%errorValues = zeros(length(range),4);
piv = 0;
for i=range
    time = dt*i;
    piv = piv+1;

    nexttile(1)
    %scatter(x_twin(:,i), y_twin(:,i), [], N_twin(:,i), 'filled');
    h = scatter(x_twin(:,i), y_twin(:,i), 2, E_twin(:,i), 'filled');
    xlim([0 20]), ylim([0 15])
    colorbar%, clim([0 3])
    clim([0 3])
    title("Individuals")

    if plotDensField
        nexttile(2)
        dFieldT = reshape(dens_twin(:,i), dims(1), dims(2));
        pcolor(dFieldT'), shading flat, colorbar
        clim([0 50])
        title('Density field')
    end


    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v);
