
direc = "D:/work/ibm-enkf/test1/"
runs3 = ["d_run_2026_4_resample_", "run_2026_4_resample_", "perturb_7_"]
prefix = strcat(direc, runs3(2));

x_e = dlmread(prefix+"e1X.csv");
y_e = dlmread(prefix+"e1Y.csv");
u_e = dlmread(prefix+"e1VX.csv");
v_e = dlmread(prefix+"e1VY.csv");

steps = -10+[99 100 101 104];

figure, tiledlayout(length(steps),3,'TileSpacing','compact')

for stepI=1:length(steps)

    
    stepVec = [x_e(:,steps(stepI))-x_e(:,steps(stepI)-1) y_e(:,steps(stepI))-y_e(:,steps(stepI)-1)];
    stepLen = sqrt(stepVec(:,1).^2 + stepVec(:,2).^2);
    spds = sqrt(u_e(:,steps(stepI)).^2 + v_e(:,steps(stepI)).^2);
    
    nexttile, histogram(spds), title('Speeds')
    nexttile, histogram(stepLen), title('Step lengths')
    nexttile, scatter(spds,stepLen, 0.4), xlabel('Speeds'), ylabel('Step lengths')


    % for i=1:50
    %     plot(x_e(i,(steps(stepI)-1):steps(stepI)),y_e(i,(steps(stepI)-1):steps(stepI)),'. -')
    %     hold on
    % end

end