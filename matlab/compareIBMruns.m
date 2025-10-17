direc = "C:/temp/";


%runs = ["d_r9_resample_", "r9_", "r9_resample_"];
runs = ["d_r14_resample_", "r14_", "r14_resample_"];
runNames = ["Dry run", "Assim", "Assim resample"];

dims_and_int = dlmread(direc+runs(1)+"fieldDims.csv");
dims = dims_and_int(1:2);
assimInt = dims_and_int(3);

for i=1:length(runs)
    prefix = strcat(direc, runs(i))

    % Read twin and ensemble densities:
    dens_twin = dlmread(strcat(prefix, "twinDens.csv"));
    dens_e = dlmread(strcat(prefix, "eDens.csv"));
    %E_twin = dlmread(strcat(prefix, "twinE.csv"));
    %E_1 = dlmread(strcat(prefix, "e1E.csv"));
    E_twin = dlmread(prefix+"twinEnergy.csv");
    E_e = dlmread(prefix+"eEnergy.csv");
    enkfField = dlmread(prefix+"enkfField.csv");

    % On first iteration, initialize arrays:
    if i==1
        rmsDens = zeros(size(dens_twin,2), length(runs))
        rmsE = zeros(size(dens_twin,2), length(runs))
        rmsEnkf = zeros(size(dens_twin,2), length(runs))
        
    end

    for j=1:size(dens_twin,2)
        j
        devi = dens_twin(:,j) - dens_e(:,j);
        rmsDens(j,i) = rms(devi);
        weightedE_twin = E_twin(:,j).*dens_twin(:,j);
        weightedE_e = E_e(:,j).*dens_e(:,j);
        %devi = E_twin(:,j) - E_e(:,j);
        devi = weightedE_twin - weightedE_e;
        rmsE(j,i) = rms(devi);

        devi = dens_e(:,j) - enkfField(:,j);
        rmsEnkf(j,i) = rms(devi);
    end
end
%%
figure, tiledlayout(2,3, "TileSpacing","compact");
nexttile, plot(rmsDens), legend(runNames), title('Density RMS'), grid on
nexttile, plot(rmsE), legend(runNames), title('Energy RMS'), grid on
% Plot rmsEnkf only at assimilation times (comparison is only valid at
% those times):
nexttile, plot(rmsEnkf(assimInt:assimInt:end,:)), legend(runNames), title('Update - EnKF field RMS'), grid on
%nexttile, plot(rmsEnkf(:,:)), legend(runNames), title('Update - EnKF field RMS'), grid on
nexttile, bar(mean(rmsDens,1)), title('Mean density RMS'), xticklabels(runNames)
nexttile, bar(mean(rmsE,1)), title('Mean energy RMS'), xticklabels(runNames)
