thresh = 0.1;
    
otdir = "D:/work/ibm-enkf/otexamples/";

ttvals = 2:2:84;
figure, tiledlayout('flow','TileSpacing','compact');
nexttile;
valsAt20 = [];
for tt=ttvals

    file = ["OT_example_t_"+string(tt)+".csv"]
    ot = dlmread(otdir+file);
    rowSums = sum(ot,2);
    sortedOt = zeros(size(ot));
    % Each row represents a from-state
    for rw = 1:size(ot,1)
        toVal = ot(rw,:);
        sortedOt(rw,:) = sort(toVal,'descend');
    end
    sortedOt = sortedOt(rowSums>thresh,:);
    for rw = 1:size(sortedOt,1)
        sortedOt(rw,:) = sortedOt(rw,:)/sortedOt(rw,1);
    end
    
    % toPlot = [4 8 10 20 30];
    % figure
    % for i=1:length(toPlot)
    %     plot((sortedOt(:,toPlot(i))))
    %     hold on
    %     grid on
    % end
    % legend(string(toPlot))
    
    semilogy(mean(sortedOt(:,1:25),1))

    hold on

    valsAt20 = [valsAt20; mean(sortedOt(:,25))];
end
grid on, legend(string(ttvals))
xlabel('Element number'), ylabel('Mean value')
nexttile, plot(ttvals, valsAt20), grid on
xlabel('Time'), ylabel('Mean val for element 20')