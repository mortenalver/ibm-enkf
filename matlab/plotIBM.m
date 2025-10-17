
%prefix = 'C:/temp/r9_';
prefix = 'C:/temp/r14_resample_';
%prefix = 'C:/temp/d_r12_resample_';

animate = 1;
plotInd = 25;

dims = dlmread([prefix 'fieldDims.csv']);

% Read twin values:

dt = 2*0.1;

x_twin = dlmread([prefix 'twinX.csv']);
y_twin = dlmread([prefix 'twinY.csv']);
E_twin = dlmread([prefix 'twinE.csv']);
N_twin = dlmread([prefix 'twinN.csv']);
U_twin = dlmread([prefix 'twinU.csv']);
V_twin = dlmread([prefix 'twinV.csv']);
dens_twin = dlmread([prefix 'twinDens.csv']);
energy_twin = dlmread([prefix 'twinEnergy.csv']);
Xfld_twin = dlmread([prefix 'twinXfld.csv']);

enkfField = dlmread([prefix 'enkfField.csv']);
% Read ensemble values:
x_1 = dlmread([prefix 'e1X.csv']);
y_1 = dlmread([prefix 'e1Y.csv']);
E_1 = dlmread([prefix 'e1E.csv']);
N_1 = dlmread([prefix 'e1N.csv']);
dens_e = dlmread([prefix 'eDens.csv']);
energy_e = dlmread([prefix 'eEnergy.csv']);

%%

rmsValues = zeros(size(x_twin,2),1);
v = VideoWriter('anim.avi');
v.FrameRate = 8; % Default 30
open(v);

figure('Renderer', 'painters', 'Position', [10 50 1400 800])
tiledlayout(2,5)

if animate > 0
    range = 1:size(x_twin,2)
else
    range = plotInd;
end
for i=range
    time = dt*i;

    nexttile(1)
    %scatter(x_twin(:,i), y_twin(:,i), [], N_twin(:,i), 'filled');
    scatter(x_twin(:,i), y_twin(:,i), [], E_twin(:,i), 'filled');
    xlim([0 20]), ylim([0 15])
    colorbar%, clim([0 3])
    clim([0 3])
    title("Individuals (twin) t="+string(time))

    nexttile(6)
    %scatter(x_1(:,i), y_1(:,i), [], N_1(:,i), 'filled');
    scatter(x_1(:,i), y_1(:,i), [], E_1(:,i), 'filled');
    xlim([0 20]), ylim([0 15])
    colorbar%, clim([0 3])
    clim([0 3])
    title('Individuals (ensemble member 1)')
    
    % nexttile(2)
    % dFieldT = reshape(Xfld_twin(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([0 2])
    % title('Food field (twin)')

    nexttile(2)
    dFieldT = reshape(dens_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 50])
    title('Density field (twin)')

    nexttile(3)
    dFieldE = reshape(dens_e(:,i), dims(1), dims(2));
    pcolor(dFieldE'), shading flat, colorbar
    clim([0 50])
    title('Density field (ensemble)')

    nexttile(4)
    diffField = dFieldT-dFieldE;
    rmsValues(i) = rms(diffField(:));
    pcolor(diffField'), shading flat, colorbar,
    title(['Deviation RMS: ' num2str(rmsValues(i),3)])
    clim([-20 20])

    % nexttile(5)
    % dFieldT = reshape(U_twin(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([-2 2])
    % title('Average U speed (twin)')
    nexttile(5)
    dFieldT = reshape(enkfField(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 50])
    title('EnKF field')

    nexttile(10)
    dFieldT = reshape(V_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([-2 2])
    title('Average V speed (twin)')

    nexttile(7)
    dFieldT = reshape(energy_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([1 3])
    title('Average energy level (twin)')

    nexttile(8)
    dFieldE = reshape(energy_e(:,i), dims(1), dims(2));
    pcolor(dFieldE'), shading flat, colorbar
    clim([1 3])
    title('Average energy level (ensemble)')

    nexttile(9)
    diffField = dFieldT-dFieldE;
    rmsv = rms(diffField(:));
    pcolor(diffField'), shading flat, colorbar,
    title(['Deviation RMS: ' num2str(rmsv,3)])
    clim([-3 3])




    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v);

figure, plot(rmsValues)