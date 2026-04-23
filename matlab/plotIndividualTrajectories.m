
%direc = "C:/temp/";
%direc = "D:/work/ibm-enkf/test2/"
direc = "E:/work/ibm-enkf/r1run_123/";
%runs = ["test5000_resample_", "test5000_"];
%runs3 = ["d_perturb_7_resample_", "perturb_7_resample_", "perturb_7_"]
%runs3 = ["d_run_2026_4_resample_", "run_2026_4_resample_", "run_2026_5_fuzzy_no_uv_upd__"]
%runs3 = ["d_loc6.5_uv_upd__resample_", "noloc_uv_upd__resample_","loc6.5_uv_upd__"];
runs3 = ["d_r2test123__resample_", "r2test123__resample_", "r2test123__"];%, "r1run2__"];%"twintest2__resample_"];%, "twintest2__resample_"];

%runs3 = ["d_indi2500_12_resample_", "indi2500_12_resample_", "indi2500_12_"];
%runs = ["test11_resample_"];
%runnames = ["Resampling"];
% Read twin and ensemble densities:
prefix = strcat(direc, runs3(1))
x_twin = dlmread(prefix+"twinX.csv");
y_twin = dlmread(prefix+"twinY.csv");
E_twin = dlmread(prefix+"twinE.csv");
N_twin = dlmread(prefix+"twinN.csv");

dt = 1*0.2;
ttt = dt*(1:size(x_twin,2));


%%
runs = runs3(1:3);
runnames = ["Free-run", "RS", "OT"];

nToPlot = 3;
idx = randperm(size(x_twin,1), nToPlot);
iFrom = 1;
iTo = 250;


lims = [0 20; 0 15];
e_lim = [0 3.5];
colo = colororder;
gcol = [0.7 0.7 0.7];
gcol2 = [0.3 0.3 0.3];

figure('Renderer', 'painters', 'Position', [50 50 1300 650]);%350]);
tiledlayout(2,1+length(runs),'TileSpacing','compact')

for i=1:length(idx)  
    nexttile(1)
    hold on, plot(x_twin(idx(i),iFrom:iTo), y_twin(idx(i),iFrom:iTo),'Color',colo(i,:),...
        'LineStyle', '-','Marker','.','MarkerSize',4)
    plot(x_twin(idx(i),iFrom), y_twin(idx(i),iFrom),'.','MarkerSize',15,'Color',colo(i,:))
    axis image, grid on, xlim(lims(1,:)), ylim(lims(2,:))
    title('Position (Twin)')
    xlabel('x position'), ylabel('y position'), box on

    nexttile(length(runs)+1+1)
    hold on, plot(ttt(iFrom:iTo)-ttt(iFrom), E_twin(idx(i),iFrom:iTo),'Color',colo(i,:),...
        'LineStyle', '-','Marker','.','MarkerSize',4)
    grid on, xlabel('Time'), box on, ylim(e_lim)
    title('Energy (Twin)')
end

for modI = 1:length(runs)
    prefix = strcat(direc, runs(modI));
    x_e = dlmread(prefix+"e1X.csv");
    y_e = dlmread(prefix+"e1Y.csv");
    E_e = dlmread(prefix+"e1E.csv");
    N_e = dlmread(prefix+"e1N.csv");
    for i=1:length(idx)  

        % Find closest particle to idx(i):
        initPos = [x_twin(idx(i),iFrom) y_twin(idx(i),iFrom)];
        initPosE = [x_e(:,iFrom) y_e(:,iFrom)];
        distSq = (initPos(1)-initPosE(:,1)).^2 + (initPos(2)-initPosE(:,2)).^2;
        [minDist, idx2] = mink(distSq, 1);

        nexttile(1+modI)
        hold on, plot(x_e(idx2,iFrom:iTo), y_e(idx2,iFrom:iTo),'Color',colo(i,:),...
            'LineStyle', ':','Marker','.','MarkerSize',4)
        plot(x_e(idx2,iFrom), y_e(idx2,iFrom),'.','MarkerSize',15,'Color',colo(i,:))
        axis image, grid on, xlim(lims(1,:)), ylim(lims(2,:))
        title("Position ("+runnames(modI)+")");
        xlabel('x position'), ylabel('y position'), box on

        nexttile(length(runs)+2+modI), grid on
        hold on, plot(ttt(iFrom:iTo)-ttt(iFrom), E_e(idx2,iFrom:iTo),'Color',colo(i,:),...
        'LineStyle', '-','Marker','.','MarkerSize',4)
        grid on, xlabel('Time'), ylim(e_lim)
        title("Energy ("+runnames(modI)+")"), box on
        
    end
end

exportgraphics(gcf, 'ind_traj.eps')

%%
% Analyse step lengths
runs = runs3([3]);
runnames = ["Free-run", "OT"];
binEdges = 0:0.01:3; %0:25;
binC = binEdges(1:end-1) + 0.5*diff(binEdges);
counts = zeros(length(binC), length(runs));
for modI = 1:length(runs)
    prefix = strcat(direc, runs(modI));
    X_e = dlmread(prefix+"e1X.csv");
    Y_e = dlmread(prefix+"e1Y.csv");

    stepLen = zeros((size(X_e,2)-1)*size(X_e,1),1);
    piv = 0;
    for i=2:size(X_e,2)
        if mod(i,assimInt) == 0%assimInt-1
            for j=1:size(X_e,1)
                piv = piv+1;
                stepVec = [X_e(j,i)-X_e(j,i-1) Y_e(j,i)-Y_e(j,i-1)];
                stepLen(piv) = sqrt(stepVec(1)^2 + stepVec(2)^2);
            end
        end
    end
    stepLen = stepLen(1:piv);
    hc = histcounts(stepLen, binEdges);
    counts(:,modI) = hc'/sum(hc);
end
figure, plot(binC, counts,'. -'), grid on
xlabel('Step size (m)'), ylabel('Relative frequency')
if length(runs) > 1
    legend(runnames)
end
exportgraphics(gcf, 'OT_stepsizes.eps')

%%
runs = runs3(2:3);
runnames = ["RS", "OT"];
%varToAnalyze = "E"; titleLabel = "Energy";
%varToAnalyze = "X"; titleLabel = "X position";
varsToAnalyze = ["X", "Y", "E"]; 
titleLabels = ["X position", "Y position", "Energy"];

perio = figure('Renderer', 'painters', 'Position', [50 50 1100 450]);%350]);
tiledlayout(1,length(varsToAnalyze));

for varI = 1:length(varsToAnalyze)
    
    varToAnalyze = varsToAnalyze(varI);
    titleLabel = titleLabels(varI);

    % figure, tiledlayout(2,1+length(runs),"TileSpacing","compact")
    perioSum = [];
    E_twin = dlmread(prefix+"twin"+varToAnalyze+".csv");
    for i=1:size(E_twin,1)
    
        t1 = E_twin(i,:);
        [psd_t1, w1] = periodogram(t1, hann(length(t1)));
        if isempty(perioSum)
            perioSum = psd_t1;
        else
            perioSum = perioSum + psd_t1;
        end

        % if rand(1) > 0.975
        %     nexttile(1), plot(t1), hold on
        %     nexttile(length(runs)+1+1), semilogy(w1, psd_t1), hold on
        % end
    end
    perioMean = perioSum/size(E_twin,1);
    % nexttile(length(runs)+1+1), plot(w1, perioMean,'k', 'LineWidth',2)
    
    perioMeanRuns = zeros(length(perioMean), length(runs));
    
    for modI = 1:length(runs)
        perioSum = [];
        prefix = strcat(direc, runs(modI));
        E_e = dlmread(prefix+"e1"+varToAnalyze+".csv");
        for i=1:size(E_e,1)
            t1 = E_e(i,:);
            [psd_t1, w1] = periodogram(t1, hann(length(t1)));
            if isempty(perioSum)
                perioSum = psd_t1;
            else
                perioSum = perioSum + psd_t1;
            end

            % if rand(1) > 0.975
            %     nexttile(1+modI), plot(t1), hold on
            %     nexttile(1+length(runs)+1+modI), semilogy(w1, psd_t1), hold on
            % end
        end
        perioMeanRuns(:,modI) = perioSum'/size(E_e,1);
        % plot(w1, perioMeanRuns(:,modI),'k', 'LineWidth',2)
    
    end
    
    % Normalized frequencies in w1 are in rad/sample.
    % To convert to frequency in (time unit)^-1:  w1 / (2 pi) / dt
    freq = w1/(2*pi)/dt;
    figure(perio), nexttile(varI)
    semilogy(freq,[perioMean perioMeanRuns])
    xlim([0 freq(end)]);
    xlabel('Frequency (d^{-1})')
    grid on; 
    if varI==1
        lgd = legend(["Twin" runnames])
        lgd.Layout.Tile = 'East'; 
    end
    title(titleLabel+" periodogram")
end

exportgraphics(gcf, 'frequencies.eps')

%%
% Find trajectories starting at almost the same position, and compare:
runs = runs3([1 3]);
runnames = ["Free-run", "OT"];
%compareRun = 1;
    
nToPlot = 3;
nPartsPer = 3;

co = colororder;

figure('Renderer', 'painters', 'Position', [50 50 1000 550]);
tiledlayout(2,nToPlot,'TileSpacing','tight')
idxAll = randperm(size(x_twin,1), nToPlot);
for compareRun=1:2
    prefix = strcat(direc, runs(compareRun));
    X_e = dlmread(prefix+"e1X.csv");
    Y_e = dlmread(prefix+"e1Y.csv");
    for iii=1:nToPlot
        % Choose one random trajectory:
        idx = idxAll(iii);
        
        fromI = 1;
        toI = 80;%55;
        initPos = [x_twin(idx,fromI) y_twin(idx,fromI)];
        
        initPosE = [X_e(:,fromI) Y_e(:,fromI)];
        distSq = (initPos(1)-initPosE(:,1)).^2 + (initPos(2)-initPosE(:,2)).^2;
        [minDist, idx2] = mink(distSq, 40);
        
        nexttile
        plot(x_twin(idx,fromI:toI), y_twin(idx,fromI:toI),'k -','LineWidth',1.5,...
            'DisplayName',"Twin")
        hold on
        %plot(x_twin(idx,fromI), y_twin(idx,fromI),'k .', 'MarkerSize', 15,...
        %    'HandleVisibility','off')
        plot(x_twin(idx,toI), y_twin(idx,toI),'k .', 'MarkerSize', 15,...
            'HandleVisibility','off')
        % Plot full trajectories for the first few:
        for i=1:nPartsPer  
            cc = [1 1 1]*(0.2 + 0.6*i/length(idx2));
            hl = plot(X_e(idx2(i),fromI:toI), Y_e(idx2(i),fromI:toI),'-', ...
                'Color',co(i,:),'DisplayName',"OT ind "+num2str(i));
            plot(X_e(idx2(i),toI), Y_e(idx2(i),toI),'.', 'MarkerSize', 15,...
                'Color',get(hl, 'Color'),'HandleVisibility','off')
        end
        % Plot end points for all:
        %plot(X_e(idx2,toI), Y_e(idx2,toI),'.','Color',0.6*[1 1 1]);
        axis image
        xlim([0 20])
        ylim([0 15])
        % xl = xlim; yl = ylim;
        % if (xl(2)-xl(1)) < (yl(2)-yl(1))
        %     mx = mean(xl);
        %     scy = 0.5*(yl(2)-yl(1));
        %     xlim([mx-scy mx+scy])
        % else
        %     my = mean(yl);
        %     scx = 0.5*(xl(2)-xl(1));
        %     ylim([my-scx my+scx])
        % end
        grid on%, xlim(lims(1,:)), ylim(lims(2,:))
        xlabel('x position'), ylabel('y position')
        title(runnames(compareRun))

        if compareRun==1 & iii==1
            lgd = legend();
            lgd.Layout.Tile = 'East'; 
        end
    end
end


% xl = [1 11]
% yl = [2 10]
% nexttile(1), xlim(xl), ylim(yl), nexttile(4), xlim(xl), ylim(yl), 
% xl = [1 11]
% yl = [6 14]
% nexttile(2), xlim(xl), ylim(yl), nexttile(5), xlim(xl), ylim(yl), 
% xl = [1 11]
% yl = [8 16]
% nexttile(3), xlim(xl), ylim(yl), nexttile(6), xlim(xl), ylim(yl), 

exportgraphics(gcf, 'ind_traj_pos0.eps')


%%
% Go through n random individuals from twin. For each, choose the closest
% individual in other simulations, and compare how far away from the twin
% they end after a certain time. Make statistics over all tested
% individuals.

runs = runs3([1 3]);
runnames = ["Free-run",  "OT"];
nInds = 5000;
fromIAll = [1 30 60]*(1/dt);
toIAll = [30 60 90]*(1/dt);%55;
titles = ["0-30 d", "30-60 d", "60-90 d"];

idxAll = randperm(size(x_twin,1), nInds)

figure('Renderer', 'painters', 'Position', [50 50 950 350]);
tiledlayout(1, 3,'TileSpacing','compact')


for interva = 1:length(fromIAll)

    fromI = fromIAll(interva);
    toI = toIAll(interva);
    endDists = zeros(nInds, length(runs));
    
    for compareRun=1:length(runs)
        prefix = strcat(direc, runs(compareRun));
        X_e = dlmread(prefix+"e1X.csv");
        Y_e = dlmread(prefix+"e1Y.csv");
        for iii=1:nInds
            % Choose one trajectory:
            idx = idxAll(iii);
            % Find init position:
            initPos = [x_twin(idx,fromI) y_twin(idx,fromI)];
            % Find the closest init position from this run:
            initPosE = [X_e(:,fromI) Y_e(:,fromI)];
            distSq = (initPos(1)-initPosE(:,1)).^2 + (initPos(2)-initPosE(:,2)).^2;
            [minDist, idx2] = min(distSq);
            if sqrt(minDist) < 0.5 % Only do this comparison if the init positions are reasonably close:
                endDists(iii,compareRun) = sqrt((X_e(idx2, toI)-x_twin(idx,toI)).^2 + (Y_e(idx2, toI)-y_twin(idx,toI)).^2) - sqrt(minDist);
            else
                endDists(iii,compareRun) = NaN;
            end
        end
    
    end
    
    edges = 0:0.25:15;
    hcAll = zeros(length(edges)-1, length(runs));
    for i=1:length(runs)
        [hc] = histcounts(endDists(:,i), edges);
        hcAll(:,i) = hc;
    end
    
    %figure%,nexttile, boxchart(endDists)
    nexttile, 
    plot(edges(1:end-1) + 0.5*diff(edges), hcAll./sum(hcAll,1),'. -') %normaliser per kolonne
    grid on, 
    if interva==1
        lgd = legend(runnames)
        lgd.Layout.Tile = 'East'; 
    end
    xlabel('Distance deviation (m)'), ylabel('Relative frequency')
    title(titles(interva))
end

exportgraphics(gcf, 'ind_traj_dev.eps')