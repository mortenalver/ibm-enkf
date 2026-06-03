
%prefix = 'C:/temp/d_salm2_';
prefix = 'C:/temp/d_test_';

animate = 1;
plotInd = 25;

dims = dlmread([prefix 'fieldDims.csv']);

% Read twin values:

dt = dims(4);
nStatesPerInd = dims(5);


speedup = 2;

% Read twin values:
twinStates = dlmread([prefix 'twin_states.csv']);
nInd = size(twinStates,1)/nStatesPerInd;

U_twin = dlmread([prefix 'twinU.csv']);
V_twin = dlmread([prefix 'twinV.csv']);
dens_twin = dlmread([prefix 'twinDens.csv']);
energy_twin = dlmread([prefix 'twinEnergy.csv']);
Xfld_twin = dlmread([prefix 'twinXfld.csv']);

% Read ensemble values:
e1States = dlmread([prefix 'e1_states.csv']);
dens_e = dlmread([prefix 'eDens.csv']);
densStd_e = dlmread([prefix 'eDensStd.csv']);
energy_e = dlmread([prefix 'eEnergy.csv']);
Xfld_e = dlmread([prefix 'eXfld.csv']);

enkfField = dlmread([prefix 'enkfField.csv']);

xMax = 8;
yMax = 8;

%%


v = VideoWriter('anim','MPEG-4');
v.FrameRate = 4; % Default 30
open(v);

figure('Renderer', 'painters', 'Position', [10 50 1100 800])
ncol = 4;%3;
tiledlayout(2,ncol)

if animate > 0
    range = 1:speedup:size(twinStates,2)
else
    range = plotInd;
end

rmsValues = zeros(length(range),4);
%errorValues = zeros(length(range),4);
piv = 0;
for i=range
    time = dt*i;
    piv = piv+1;

    twinSt = reshape(twinStates(:,i),[nStatesPerInd nInd]);
    x_twin = twinSt(1,:);
    y_twin = twinSt(2,:);
    vx_twin = twinSt(3,:);
    vy_twin = twinSt(4,:);
    E_twin = twinSt(6,:);

    e1St = reshape(e1States(:,i),[nStatesPerInd nInd]);
    x_1 = e1St(1,:);
    y_1 = e1St(2,:);
    vx_1 = e1St(3,:);
    vy_1 = e1St(4,:);
    E_1 = e1St(6,:);

    nexttile(1)
    hold off
    
    %scatter(x_twin, y_twin, [], N_twin(:,i), 'filled');
    h = scatter(x_twin, y_twin, 3, E_twin, 'filled');
    xlim([0 xMax]), ylim([0 yMax]), , axis image
    %colorbar%, clim([0 3])
    clim([0 0.1])
    rectangle('Position',[0 0 8 8],'Curvature',[1 1])
    title("Individuals (twin) t="+string(time))
    axis image, xlim([0 8]), ylim([0 8]);
    
    nexttile(ncol+1)
    hold off
    scatter(x_1, y_1, 3, E_1, 'filled');
    clim([0 0.1])
    %colorbar, clim([0 3])
    rectangle('Position',[0 0 8 8],'Curvature',[1 1])
    title('Individuals (ensemble member 1)')
    axis image, xlim([0 8]), ylim([0 8]);

    nexttile(4)
    dFieldT = reshape(Xfld_e(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 1])
    title('Food field (ensemble)')

    nexttile(ncol+4)
    dFieldT = reshape(Xfld_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 1])
    title('Food field (twin)')

    % nexttile(ncol+4)
    % dFieldT = reshape(densStd_e(:,i), dims(1), dims(2));
    % pcolor(dFieldT'), shading flat, colorbar
    % clim([0 10])
    % title('Density standard dev. (ensemble)')

    nexttile(2)
    dFieldT = reshape(dens_twin(:,i), dims(1), dims(2));
    pcolor(dFieldT'), shading flat, colorbar
    clim([0 6])
    title('Density field (twin)')

    nexttile(ncol+2)
    dFieldE = reshape(dens_e(:,i), dims(1), dims(2));
    pcolor(dFieldE'), shading flat, colorbar
    clim([0 6]); %clim([0 50])
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