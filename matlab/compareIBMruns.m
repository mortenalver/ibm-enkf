%direc = "C:/temp/";
%direc = "D:/work/ibm-enkf/test1/"
%direc = "D:/work/ibm-enkf/test2/"
%direc = ["D:/work/ibm-enkf/r1run/", "D:/work/ibm-enkf/r1run/", "D:/work/ibm-enkf/r1run/"];
direc = ["C:/temp/", "C:/temp/", "C:/temp/", "C:/temp/", "C:/temp/"];
%direc = ["D:/work/ibm-enkf/r1run/", "D:/work/ibm-enkf/r1run/", "D:/work/ibm-enkf/r1run/"];
%runs = ["d_test11_resample_", "test12_resample_", "test12_"];

% runs = ["d_test2500_resample_", "test2500_resample_", "test2500_3_resample_", "test2500_4_resample_", "test2500_5_resample_"]%, "test5000_"];
% runNames = ["Free-run", "Resample", "Resample 3", "Resample 4", "Resample 5"]%, "Optimal transport"];

%runs = ["d_perturb_5_resample_", "perturb_7_resample_", "perturb_assimint6_resample_", "perturb_7_", "perturb_assimint6_"];
%runs = ["d_run_2026_4_resample_", "run_2026_4_resample_", "run_2026_4_", "run_2026_4_no_uv_upd__", "run_2026_5_small_uv_upd__"];
%runs = ["d_run_2026_6_no_uv_upd__resample_", "run_2026_6_no_uv_upd__resample_","run_2026_7_no_uv_upd__","run_2026_8_no_uv_upd__"];
%runs = ["d_higheps_1_resample_", "higheps_1_resample_", "higheps_2_"];
%runs = ["d_loc6.5_uv_upd__resample_", "loc6.5_uv_upd__resample_","loc6.5_uv_upd__"];%"loc1.5_uv_upd__resample_", "loc_1_no_uv_upd__resample_"];

%runNames = ["Free-run", "RS", "OT"];%"RS 1.5", "RS 2.5", %, "OT 2"];%, "OT 3"];%, "OT rapid"];
%runs = ["d_twintest2__resample_", "twintest2__resample_", "twintest2__"];

%runs = ["d_r1run__resample_", "r1run2__resample_";];%"twintest2__resample_"];%, "twintest2__resample_"];
%runNames = ["Free-run", "RS"];%, "RS 20"];

%runs = ["d_test1__resample_", "test1__resample_", "test_nouv__resample_", "test_nospeed__resample_"];%, ...
    %"anamorph__resample_"];
%runNames = ["Free-run", "RS", "RS no food", "RS no speed"];%, "RS GAT"];%, "RS 20"];
runs = ["d_test1__resample_", "test1__resample_", "testr__resample_"];
runNames = ["Free-run", "RS", "RS lower R"];%, "RS 20"];


markSignificances = 0;

plotEndTime = 10;%100;

dt = 1*0.2;

plotDists = 1;
plotTimes = [10 20 35];

signLevel = 0.05;


dims_and_int = dlmread(direc(1)+runs(1)+"fieldDims.csv");
dims = dims_and_int(1:2);
assimInt = dims_and_int(3);

if plotDists > 0
    figure('Renderer', 'painters', 'Position', [50 50 850 500]);
    tiledlayout(length(plotTimes), 1+length(runs),'TileSpacing','compact')
end

for i=1:length(runs)
    prefix = strcat(direc(i), runs(i))

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
        %corrDens2 = zeros(size(dens_twin,2), length(runs));
        corrE = zeros(size(dens_twin,2), length(runs));
        corrX = zeros(size(dens_twin,2), length(runs));
    end

    for j=1:size(dens_twin,2)
        
        devi = dens_twin(:,j) - dens_e(:,j);
        rmsDens(j,i) = myRms(devi,dims);
        stdDens(j,i) = mean(densStd_e(:,j));
        corrval = corrcoef([dens_twin(:,j) dens_e(:,j)]);
        if ~isnan(corrval(1,2))
            corrDens(j,i) = corrval(1,2);
        end
        % corrDens(j,i) = smoothedCorrelation(dens_twin(:,j), dens_e(:,j), dims, 0);
        
        weightedE_twin = E_twin(:,j).*dens_twin(:,j);
        weightedE_e = E_e(:,j).*dens_e(:,j);
        %devi = E_twin(:,j) - E_e(:,j);
        devi = weightedE_twin - weightedE_e;
        rmsE(j,i) = myRms(devi,dims);
        corrval = corrcoef([weightedE_twin weightedE_e]);
        if ~isnan(corrval(1,2))
            corrE(j,i) = corrval(1,2);
        end
        % corrE(j,i) = smoothedCorrelation(weightedE_twin, weightedE_e, dims, 0);
        
        devi = Xfld_twin(:,j) - Xfld_e(:,j);
        rmsX(j,i) = myRms(devi,dims);
        corrval = corrcoef([Xfld_twin(:,j) Xfld_e(:,j)]);
        if ~isnan(corrval(1,2))
            corrX(j,i) = corrval(1,2);
        end
        %corrX(j,i) = smoothedCorrelation(Xfld_twin(:,j), Xfld_e(:,j), dims, 0);

        devi = dens_e(:,j) - enkfField(:,j);
        rmsEnkfIBM(j,i) = myRms(devi,dims);
        devi = dens_twin(:,j) - enkfField(:,j);
        rmsEnkf(j,i) = myRms(devi, dims);
        
        if plotDists > 0
            clims = [0 20];%[0 3];%[0 100];
            ppp = find(j*dt == plotTimes);
            if numel(ppp)>0
                j
                if i==1
                    % Plot twin:
                    nexttile((ppp-1)*(1+length(runs))+1)
                    dField = reshape(dens_twin(:,j), dims(1), dims(2));
                    %dField = reshape(Xfld_twin(:,j), dims(1), dims(2));
                    pcolor(dField'), shading flat, axis image
                    if ppp==1
                        cbr = colorbar, 
                        cbr.Layout.Tile = 'East'; 
                    end
                    clim(clims)
                    %scatter(x_twin(:,i), y_twin(:,j), 1, Ei_twin(:,i), 'filled');
                    %xlim([0 20]), ylim([0 15]), colorbar, clim([0 3])
                    title("Twin (t="+string(plotTimes(ppp))+")")
                    axis off, box on
                end
                nexttile((ppp-1)*(1+length(runs))+1+i)
                dField = reshape(dens_e(:,j), dims(1), dims(2));
                %dField = reshape(Xfld_e(:,j), dims(1), dims(2));
                pcolor(dField'), shading flat, axis image
                clim(clims)
                %scatter(x_1(:,i), y_1(:,i), 1, Ei1(:,j), 'filled');
                %xlim([0 20]), ylim([0 15]), colorbar, clim([0 3])
                title(runNames(i)+" (t="+string(plotTimes(ppp))+")")
                axis off, box on
            end
        end
    end
end
exportgraphics(gcf, 'density_snaps.eps')
%%
% Set ttt in case the other matrices have grown larger:
ttt = dt*(1:size(rmsDens,1));

figure('Renderer', 'painters', 'Position', [50 50 1000 650]);
tld = tiledlayout(4,3, "TileSpacing","compact");
nexttile, plot(ttt, rmsDens), title('Density RMS'), grid on, hold on
xlabel('Time')
lgd = legend([runNames])%, lgd.Location = 'NorthEast'
lgd.Layout.Tile = 'East'; 
fixylim(gca, plotEndTime);
nexttile, plot(ttt, rmsE), title('Energy RMS'), grid on
xlabel('Time')

nexttile, plot(ttt, rmsX), title('Food field RMS'), grid on
xlabel('Time')

nexttile, bar(mean(rmsDens,1)), title('Mean density RMS'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(rmsDens,1), std(rmsDens,0,1),'Color','k','LineStyle','none','LineWidth',1);
if markSignificances
    markSignificance(h, runs, rmsDens, signLevel)
end
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
if markSignificances
    markSignificance(h, runs, rmsE, signLevel)
end

nexttile, bar(mean(rmsX,1)), title('Mean food field RMS'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(rmsX,1), std(rmsX,0,1),'Color','k','LineStyle','none','LineWidth',1);
if markSignificances
    markSignificance(h, runs, rmsX, signLevel)
end

nexttile, plot(ttt, corrDens), title('Density correlation'), grid on, hold on
xlabel('Time')

nexttile, plot(ttt, corrE), title('Energy correlation'), grid on
xlabel('Time')

nexttile, plot(ttt, corrX), title('Food field correlation'), grid on
xlabel('Time')


nexttile, bar(mean(corrDens,1)), title('Mean density correlation'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(corrDens,1), std(corrDens,0,1),'Color','k','LineStyle','none','LineWidth',1);
if markSignificances
    markSignificance(h, runs, corrDens, signLevel)
end
ylim([0 1])

nexttile, bar(mean(corrE,1)), title('Mean energy correlation'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(corrE,1), std(corrE,0,1),'Color','k','LineStyle','none','LineWidth',1);
if markSignificances
    markSignificance(h, runs, corrE, signLevel)
end
ylim([0 1])

nexttile, bar(mean(corrX,1)), title('Mean food field correlation'), xticklabels(runNames), grid on
hold on, h = errorbar(1:length(runs), mean(corrX,1), std(corrX,0,1),'Color','k','LineStyle','none','LineWidth',1);
if markSignificances
    markSignificance(h, runs, corrX, signLevel)
end
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
if markSignificances
    markSignificance(h, runs(inclR), rmsEnkf(assimInt:assimInt:end,inclR), signLevel)
end
xticklabels(runNames(inclR)), grid on
title('EnKF accuracy'), grid on
yl = ylim;

nexttile, bar(mean(rmsEnkfIBM(assimInt:assimInt:end,inclR)))
hold on, h = errorbar(1:length(inclR), mean(rmsEnkfIBM(assimInt:assimInt:end,inclR),1), std(rmsEnkfIBM(assimInt:assimInt:end,inclR),0,1),'Color','k','LineStyle','none','LineWidth',1);
if markSignificances
    markSignificance(h, runs(inclR), rmsEnkfIBM(assimInt:assimInt:end,inclR), signLevel)
end
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
            text(h.XData(i)+0.1, (bar_height), '*','FontSize',18);
        end
    end

    for i=3:length(runs)
        aov = anova(values(:,[2 i]));
        pval = aov.stats.pValue(1);
        if pval < signLevel
            bar_height = h.YData(i)*1.075;%+ 0.15* h.YPositiveDelta(i);
            text(h.XData(i)+0.15, (bar_height), ' ^o','FontSize',16);
        end
    end

end

function fixylim(gca, plotEndTime)
    yl = ylim;
    if yl(2) > plotEndTime
        ylim([0 plotEndTime]);
    end

end