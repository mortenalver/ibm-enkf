%direc = "C:/temp/";
direc = "D:/work/ibm-enkf/test1/"

%runs = ["d_test11_resample_", "test12_resample_", "test12_"];

% runs = ["d_test2500_resample_", "test2500_resample_", "test2500_3_resample_", "test2500_4_resample_", "test2500_5_resample_"]%, "test5000_"];
% runNames = ["Free-run", "Resample", "Resample 3", "Resample 4", "Resample 5"]%, "Optimal transport"];

runs = ["d_perturb_3_resample_", "perturb_3_resample_", "perturb_3_"];
%runs = ["d_test_1_resample_", "test_1_resample_", "test_1_"];%"indi2500_12_"];
runNames = ["Free-run", "RS", "OT"];


dt = 2*0.2;

plotDists = 1;
plotTimes = [10 50 100];

signLevel = 0.05;


dims_and_int = dlmread(direc+runs(1)+"fieldDims.csv");
dims = dims_and_int(1:2);
assimInt = dims_and_int(3);

if plotDists > 0
    figure, tiledlayout(length(plotTimes), 1+length(runs),'TileSpacing','compact')
end

for i=1:length(runs)
    prefix = strcat(direc, runs(i))

    % Read twin and ensemble densities:
    dens_twin = dlmread(strcat(prefix, "twinDens.csv"));
    dens_e = dlmread(strcat(prefix, "eDens.csv"));
    densStd_e = dlmread(strcat(prefix, "eDensStd.csv"));
    Xfld_twin = dlmread(strcat(prefix, "twinXfld.csv"));
    Xfld_e = dlmread(strcat(prefix, "eXfld.csv"));
    %E_twin = dlmread(strcat(prefix, "twinE.csv"));
    %E_1 = dlmread(strcat(prefix, "e1E.csv"));
    E_twin = dlmread(prefix+"twinEnergy.csv");
    E_e = dlmread(prefix+"eEnergy.csv");
    enkfField = dlmread(prefix+"enkfField.csv");
    x_twin = dlmread(prefix+"twinX.csv");
    y_twin = dlmread(prefix+"twinY.csv");
    Ei_twin = dlmread(prefix+"twinE.csv");
    x_1 = dlmread(prefix+"e1X.csv");
    y_1 = dlmread(prefix+"e1Y.csv");
    Ei1 = dlmread(prefix+"e1E.csv");

    % On first iteration, initialize arrays:
    if i==1
        ttt = dt*(1:size(dens_twin,2));
        rmsDens = zeros(size(dens_twin,2), length(runs));
        rmsE = zeros(size(dens_twin,2), length(runs));
        rmsX = zeros(size(dens_twin,2), length(runs));
        rmsEnkf = zeros(size(dens_twin,2), length(runs));
        rmsEnkfIBM = zeros(size(dens_twin,2), length(runs));
        stdDens = zeros(size(dens_twin,2), length(runs));
        corrDens = zeros(size(dens_twin,2), length(runs));
        corrE = zeros(size(dens_twin,2), length(runs));
        corrX = zeros(size(dens_twin,2), length(runs));
    end

    for j=1:size(dens_twin,2)
        
        devi = dens_twin(:,j) - dens_e(:,j);
        rmsDens(j,i) = rms(devi);
        stdDens(j,i) = mean(densStd_e(:,j));
        corrval = corrcoef([dens_twin(:,j) dens_e(:,j)]);
        if ~isnan(corrval(1,2))
            corrDens(j,i) = corrval(1,2);
        end
        weightedE_twin = E_twin(:,j).*dens_twin(:,j);
        weightedE_e = E_e(:,j).*dens_e(:,j);
        %devi = E_twin(:,j) - E_e(:,j);
        devi = weightedE_twin - weightedE_e;
        rmsE(j,i) = rms(devi);
        corrval = corrcoef([weightedE_twin weightedE_e]);
        if ~isnan(corrval(1,2))
            corrE(j,i) = corrval(1,2);
        end

        devi = Xfld_twin(:,j) - Xfld_e(:,j);
        rmsX(j,i) = rms(devi);
        corrval = corrcoef([Xfld_twin(:,j) Xfld_e(:,j)]);
        if ~isnan(corrval(1,2))
            corrX(j,i) = corrval(1,2);
        end

        devi = dens_e(:,j) - enkfField(:,j);
        rmsEnkfIBM(j,i) = rms(devi);
        devi = dens_twin(:,j) - enkfField(:,j);
        rmsEnkf(j,i) = rms(devi);
        
        if plotDists > 0
            clims = [0 100];%[0 3];%[0 100];
            ppp = find(j*dt == plotTimes);
            if numel(ppp)>0
                j
                if i==1
                    % Plot twin:
                    nexttile((ppp-1)*(1+length(runs))+1)
                    dField = reshape(dens_twin(:,j), dims(1), dims(2));
                    %dField = reshape(Xfld_twin(:,j), dims(1), dims(2));
                    pcolor(dField'), shading flat, colorbar, axis image
                    clim(clims)
                    %scatter(x_twin(:,i), y_twin(:,j), 1, Ei_twin(:,i), 'filled');
                    %xlim([0 20]), ylim([0 15]), colorbar, clim([0 3])
                    title("Twin (t="+string(plotTimes(ppp))+")")
                end
                nexttile((ppp-1)*(1+length(runs))+1+i)
                dField = reshape(dens_e(:,j), dims(1), dims(2));
                %dField = reshape(Xfld_e(:,j), dims(1), dims(2));
                pcolor(dField'), shading flat, colorbar, , axis image
                clim(clims)
                %scatter(x_1(:,i), y_1(:,i), 1, Ei1(:,j), 'filled');
                %xlim([0 20]), ylim([0 15]), colorbar, clim([0 3])
                title(runNames(i)+" (t="+string(plotTimes(ppp))+")")
            end
        end
    end
end
%%
figure('Renderer', 'painters', 'Position', [50 50 1000 650]);
tld = tiledlayout(4,3, "TileSpacing","compact");
nexttile, plot(ttt, rmsDens), title('Density RMS'), grid on, hold on
xlabel('Time')
lgd = legend([runNames]), lgd.Location = 'NorthEast'


nexttile, plot(ttt, rmsE), title('Energy RMS'), grid on
xlabel('Time')

nexttile, plot(ttt, rmsX), title('Food field RMS'), grid on
xlabel('Time')

nexttile, bar(mean(rmsDens,1)), title('Mean density RMS'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(rmsDens,1), std(rmsDens,0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs, rmsDens, signLevel)
% for i=2:length(runs)
%     aov = anova(rmsDens(:,[1 i]));
%     pval = aov.stats.pValue(1);
%     if pval < 0.05
%         bar_height = h.YData(i)+ 0.15* h.YPositiveDelta(i);
%         text(h.XData(i)+0.1, (bar_height), '*','FontSize',18);
%     end
% end
%sigstar(groups, pvals)

nexttile, bar(mean(rmsE,1)), title('Mean energy RMS'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(rmsE,1), std(rmsE,0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs, rmsE, signLevel)

nexttile, bar(mean(rmsX,1)), title('Mean food field RMS'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(rmsX,1), std(rmsX,0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs, rmsX, signLevel)

nexttile, plot(ttt, corrDens), title('Density correlation'), grid on, hold on
xlabel('Time')

nexttile, plot(ttt, corrE), title('Energy correlation'), grid on
xlabel('Time')

nexttile, plot(ttt, corrX), title('Food field correlation'), grid on
xlabel('Time')


nexttile, bar(mean(corrDens,1)), title('Mean density correlation'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(corrDens,1), std(corrDens,0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs, corrDens, signLevel)
ylim([0 1])

nexttile, bar(mean(corrE,1)), title('Mean energy correlation'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(corrE,1), std(corrE,0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs, corrE, signLevel)
ylim([0 1])

nexttile, bar(mean(corrX,1)), title('Mean food field correlation'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(corrX,1), std(corrX,0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs, corrX, signLevel)
ylim([0 1])


exportgraphics(gcf, 'run_stats.eps')

%%
inclR = 2:length(runs);
figure('Renderer', 'painters', 'Position', [50 50 800 300]);
tiledlayout(1,2, "TileSpacing","compact")
%nexttile, plot(rmsEnkfIBM(assimInt:assimInt:end,2:end)), legend(runNames,'Location','southeast')
%legend(runNames,'Location','southeast')

%nexttile, plot(rmsEnkf(:,:)), legend(runNames), title('Update - EnKF field RMS'), grid on

nexttile, bar(mean(rmsEnkf(assimInt:assimInt:end,2:end))), grid on %, legend(runNames,'Location','southeast')
hold on, h = errorbar(1:length(inclR), mean(rmsEnkf(assimInt:assimInt:end,inclR),1), std(rmsEnkf(assimInt:assimInt:end,inclR),0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs(inclR), rmsEnkf(assimInt:assimInt:end,inclR), signLevel)
xticklabels(runNames(inclR)), grid on
title('EnKF accuracy'), grid on
yl = ylim;

nexttile, bar(mean(rmsEnkfIBM(assimInt:assimInt:end,inclR)))
hold on, h = errorbar(1:length(inclR), mean(rmsEnkfIBM(assimInt:assimInt:end,inclR),1), std(rmsEnkfIBM(assimInt:assimInt:end,inclR),0,1),'Color','k','LineStyle','none','LineWidth',1);
markSignificance(h, runs(inclR), rmsEnkfIBM(assimInt:assimInt:end,inclR), signLevel)
xticklabels(runNames(inclR)), grid on
ylim(yl)
title('IBM update accuracy'), grid on

exportgraphics(gcf, 'update_stats.eps')

%%
%nexttile, plot(stdDens), title('Mean density stddev'), grid on, legend(runNames,'Location','southeast')



function markSignificance(h, runs, values, signLevel)

    for i=2:length(runs)
        aov = anova(values(:,[1 i]));
        pval = aov.stats.pValue(1);
        if pval < signLevel
            bar_height = h.YData(i)*1.075;%+ 0.15* h.YPositiveDelta(i);
            text(h.XData(i)+0.1, (bar_height), '*','FontSize',22);
        end
    end

end