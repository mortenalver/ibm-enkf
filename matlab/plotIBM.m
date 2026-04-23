
%prefix = 'D:/work/ibm-enkf/test1/perturb_5_';
%prefix = 'D:/work/ibm-enkf/test1/run_2026_3_resample_';
%prefix = 'D:/work/ibm-enkf/test2/loc6.5_uv_upd__';
%prefix = 'D:/work/ibm-enkf/r1run/twintest__resample_';
%prefix = 'C:/temp/twintest2__resample_';
prefix = 'C:/temp/d_twin123_resample_';
%prefix = 'C:/temp/d_gtwin_perturb_1_resample_';
%prefix = 'C:/temp/d_gtwin_indi2500_13_resample_';
%prefix = 'C:/temp/indi2500_10_resample_';
%prefix = 'C:/temp/nomigr11_resample_';
%prefix = 'C:/temp/d_nomigr_long1_resample_';

animate = 1;
plotInd = 25;

dims = dlmread([prefix 'fieldDims.csv']);

% Read twin values:

dt = 1*0.2;
if length(dims) > 3
    dt = dims(4);
end
speedup = 5;

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

figure('Renderer', 'painters', 'Position', [10 50 1100 800])
ncol = 4;%3;
tiledlayout(2,ncol)

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
    h = scatter(x_twin(:,i), y_twin(:,i), 1, E_twin(:,i), 'filled');
    xlim([0 20]), ylim([0 15])
    colorbar%, clim([0 3])
    clim([0 3])
    title("Individuals (twin) t="+string(time))

    nexttile(ncol+1)
    %scatter(x_1(:,i), y_1(:,i), [], N_1(:,i), 'filled');
    scatter(x_1(:,i), y_1(:,i), 1, E_1(:,i), 'filled');
    xlim([0 20]), ylim([0 15])
    colorbar%, clim([0 3])
    clim([0 3])
    title('Individuals (ensemble member 1)')
    
    nexttile(4)
    dFieldT = reshape(Xfld_e(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 2])
    title('Food field (ensemble)')

    nexttile(ncol+4)
    dFieldT = reshape(Xfld_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 2])
    title('Food field (twin)')

    % nexttile(ncol+4)
    % dFieldT = reshape(densStd_e(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([0 10])
    % title('Density standard dev. (ensemble)')

    nexttile(2)
    dFieldT = reshape(dens_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 50])
    title('Density field (twin)')

    nexttile(ncol+2)
    dFieldE = reshape(dens_e(:,i), dims(1), dims(2));
    pcolor(dFieldE'), shading flat, colorbar
    clim([0 50])
    title('Density field (ensemble)')

    nexttile(3)
    dFieldT = reshape(dens_twin(:,i), dims(1), dims(2));
    dFieldE = reshape(dens_e(:,i), dims(1), dims(2));
    diffField = dFieldT-dFieldE;
    rmsValues(piv,1) = rms(diffField(:));
    rmsValues(piv,2) = mean(abs(diffField(:)));
    sortedErr = sort(abs(diffField(:)));
    rmsValues(piv,3) = sortedErr(round(0.9*length(sortedErr)));
    corrvalues = corrcoef([dFieldT(:) dFieldE(:)]);
    rmsValues(piv,4) = corrvalues(1,2);
    pcolor(diffField'), shading flat, colorbar,
    title(['Deviation RMS: ' num2str(rmsValues(piv),3)])
    clim([-20 20])

    nexttile(ncol+3)
    histogram(diffField(:))
    title('Density deviations')

    % nexttile(5)
    % dFieldT = reshape(U_twin(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([-2 2])
    % title('Average U speed (twin)')
    % nexttile(5)
    % dFieldT = reshape(enkfField(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([0 50])
    % title('EnKF field')

    % nexttile(10)
    % dFieldT = reshape(V_twin(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([-2 2])
    % title('Average V speed (twin)')

    % nexttile(ncol+3)
    % dFieldT = reshape(energy_twin(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % %clim([1 3])
    % title('Average energy level (twin)')
    % 
    % nexttile(ncol+4)
    % dFieldE = reshape(energy_e(:,i), dims(1), dims(2));
    % pcolor(dFieldE'), shading flat, colorbar
    % clim([1 3])
    % title('Average energy level (ensemble)')

    % nexttile(9)
    % diffField = dFieldT-dFieldE;
    % rmsv = rms(diffField(:));
    % pcolor(diffField'), shading flat, colorbar,
    % title(['Deviation RMS: ' num2str(rmsv,3)])
    % clim([-3 3])




    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v);

figure
nexttile, plot(rmsValues(:,1:3)), grid on, title('Deviation metrics')
legend('RMS', 'MAE','90 percentile')
nexttile, plot(rmsValues(:,4)), grid on, title('Correlation coefficient')