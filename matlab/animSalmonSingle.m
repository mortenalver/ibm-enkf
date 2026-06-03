
prefix = 'C:/temp/d_salm2_';

animate = 1;
plotInd = 25;

dims = dlmread([prefix 'fieldDims.csv']);

% Read twin values:

dt = dims(4);
nStatesPerInd = dims(5);
speedup = 1;

enkfField = dlmread([prefix 'enkfField.csv']);
% Read ensemble member 1 values:
e1States = dlmread([prefix 'e1_states.csv']);

xMax = 8;
yMax = 8;

%%

v = VideoWriter('salmon','MPEG-4');
v.FrameRate = 4; % Default 30
open(v);

figure('Renderer', 'painters', 'Position', [10 50 1100 800])

rmsValues = zeros(length(range),4);
for i=range
    time = dt*i;

    e1St = reshape(e1States(:,i),[nStatesPerInd nInd]);
    x_1 = e1St(1,:);
    y_1 = e1St(2,:);
    vx_1 = e1St(3,:);
    vy_1 = e1St(4,:);
    E_1 = e1St(6,:);
    hold off
    for j=1:length(x_1)
        xval = [x_1(j) x_1(j)-0.5*vx_1(j)];
        yval = [y_1(j) y_1(j)-0.5*vy_1(j)];
        plot(xval, yval), hold on
    end
    h = rectangle('Position',[0 0 8 8],'Curvature',[1 1])
    set(h,"LineWidth",2)
    title("Individuals t="+string(time))
    axis image, xlim([0 8]), ylim([0 8]);
    
    

    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v);

