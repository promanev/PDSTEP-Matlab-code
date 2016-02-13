function oneStepAnalyze()
% Analyze robots that spend more time upright after one step forward vs.
% those that spend less time upright.
%
% Variables: joint angles/ang.velocity/ang.acceleration, CoM traces, joint
% kinetics (forces applied?) - this might be a comparison to EMG in humans,
% foot-surface contact forces.
% clear
% clc
% Load settings:
% if exist('DEsettings.mat','file') ~= 0
%     load('DEsettings.mat');
%     disp('...Using stored DE settings: ')
%     fdNames = fieldnames(DEsettings);
%     for i = 1:numel(fdNames)
%         disp([fdNames{i} ' = ' num2str(DEsettings.(fdNames{i}))]);
%     end
%     key = input('Are these settings OK? (y/n): ', 's');
%     if key~='y'
%         error('Initialization stopped. Change DE settings and save as DEsettings.m')
%     end
% else
%     error('NO SETTINGS FILE');
% end

% wts_size = [DEsettings.numInput, DEsettings.numHidden, DEsettings.numOutput];

fileList = dir('weights*.txt');
numFiles = size(fileList,1);
disp(['Running ' num2str(numFiles) ' files...'])

% cycle through all weight files:
for i=1:numFiles
    [jointAngs,com,jointForces,lf,rf,targets,swingFootCOMtrace, swingFootTouch]=runSim();
end




% Joint forces/angles are collected before initialized during step 1 of simulation
% therefore these data are deleted:
% jointAngs(jointAngs<-100)=0;
% jointForces = jointForces(:,2:end);

%Start plotting:
plotCOM(com,targets,swingFootCOMtrace);
plotJointAngles(jointAngs);
disp('THE END')
end



%% Misc functions:
function [jointAngs, com, jointForces, leftFootForce, rightFootForce,...
    targets, swingFootCOMtrace, swingFootTouch] = runSim(fileNameArray,execName)
% Function runs the simulation and reads output files:
if nargin<2
    execName = 'PDSTEP_analyze';
end

if nargin<1
    fileNameArray = {'jointAngs', 'com', 'jointForces',...
        'leftFootForce', 'rightFootForce','targets','swingFootCOMtrace',...
        'swingFootTouch'};
    % how many columns to extract from .txt files:
    fileDataSize = [12, 3, 12, 3, 3, 3, 3, 1];                    
end

% check that all files are present:
if exist('GLUT32.DLL','file')==0
    warning('!!! MISSING GLUT32.DLL !!!')
    input('Copy necessary file and press any key','s');
end

if exist([execName '.exe'],'file')==0
    warning(['!!! MISSING ' execName '.exe !!!'])
    input('Copy necessary file and press any key','s');
end

if exist([execName '.ilk'],'file')==0
    warning(['!!! MISSING ' execName '.ilk !!!'])
    input('Copy necessary file and press any key','s');
end

if exist([execName '.pdb'],'file')==0
    warning(['!!! MISSING ' execName '.pdb !!!'])
    input('Copy necessary file and press any key','s');
end

system([execName '.exe']);

% Making sure there are no un-closed files:
if ~isempty(fopen('all'))
    fclose('all');
end

for i=1:length(fileNameArray)
    if exist([fileNameArray{i} '.txt'],'file')~=0
        fid = fopen([fileNameArray{i} '.txt'],'r');
        varName = genvarname(fileNameArray{i});
        eval([varName '= fscanf(fid,''%f'',[fileDataSize(i),Inf]);']);
        fclose(fid);
        delete([fileNameArray{i} '.txt']);
    else
        warning([fileNameArray{i} '.txt is missing. Proceeding to the next file...']);
    end
end

end

function plotCOM(f,t,swingFootCOMtrace)
% mirror x-axis values, so the graph is the same as looking at robot's back
% as in simulation: 
f(1,:) = -1*f(1,:);
t(1,:) = -1*t(1,:);
swingFootCOMtrace(1,:) = -1*swingFootCOMtrace(1,:);

XMIN = -3.5;
XMAX = 3.5;
YMIN = -0.5;
YMAX = 3.5;
ZMIN = -0.5;
ZMAX = 3;
footLength = 0.812908;
footWidth = 0.541939;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)

plot(f(1,:),f(2,:)) % COM trace in 2D
hold on
plot(f(1,1),f(2,1),'ro') % the initial position
plot(f(1,end),f(2,end),'ko') % the end position
plot(t(1,1), t(2,1), 'gx') % left target
plot(t(1,2), t(2,2), 'bx') % right target
plot(swingFootCOMtrace(1,:),swingFootCOMtrace(2,:),'g')
plot(swingFootCOMtrace(1,1),swingFootCOMtrace(2,1),'ro') % the initial position
plot(swingFootCOMtrace(1,end),swingFootCOMtrace(2,end),'ko') % the end position
title('COM path in XY plane (BACK)')
xlabel('X-axis (LEFT \Rightarrow RIGHT)')
ylabel('Y-axis (DOWN \Rightarrow UP)')
grid on
legend('Whole-body COM', 'Start', 'End','Left Target', 'Right Target', 'Left Foot COM', 'Location', 'Best')
axis([XMIN XMAX YMIN YMAX])
hold off

subplot(1,3,2)

plot(f(1,:),f(3,:))
hold on
plot(f(1,1),f(3,1),'ro')
plot(f(1,end),f(3,end),'ko')
plot(t(1,1), t(3,1), 'gx')
plot(t(1,2), t(3,2), 'bx')
plot(swingFootCOMtrace(1,:),swingFootCOMtrace(3,:),'g')
plot(swingFootCOMtrace(1,1),swingFootCOMtrace(3,1),'ro') % the initial position
plot(swingFootCOMtrace(1,end),swingFootCOMtrace(3,end),'ko') % the end position
title('COM path in XZ plane (BIRD''S-EYE VIEW)')
xlabel('X-axis (LEFT \Rightarrow RIGHT)')
ylabel('Z-axis (BACKWARD \Rightarrow FORWARD)')
% Draw base of support:
%back line
line([(-0.3609-footWidth/2) (0.3609+footWidth/2)],[-footLength/2 -footLength/2])
% front line
line([(-0.3609-footWidth/2) (0.3609+footWidth/2)],[footLength/2 footLength/2])
% left line
line([(-0.3609-footWidth/2) (-0.3609-footWidth/2)],[-footLength/2 footLength/2])
%right line
line([(0.3609+footWidth/2) (0.3609+footWidth/2)],[-footLength/2 footLength/2])
grid on
legend('Whole-body COM', 'Start', 'End','Left Target', 'Right Target', 'Left Foot COM', 'Location', 'Best')
axis([XMIN XMAX ZMIN ZMAX])
hold off

subplot(1,3,3)

plot(f(3,:),f(2,:))
hold on
plot(f(3,1),f(2,1),'ro')
plot(f(3,end),f(2,end),'ko')
plot(t(3,1), t(2,1), 'gx')
plot(t(3,2), t(2,2), 'bx')
plot(swingFootCOMtrace(3,:),swingFootCOMtrace(2,:),'g')
plot(swingFootCOMtrace(3,1),swingFootCOMtrace(2,1),'ro') % the initial position
plot(swingFootCOMtrace(3,end),swingFootCOMtrace(2,end),'ko') % the end position
title('COM path in ZY plane (RIGHT SIDE)')
xlabel('Z-axis (BACKWARD \Rightarrow FORWARD)')
ylabel('Y-axis (DOWN \Rightarrow UP)')
grid on
legend('Whole-body COM', 'Start', 'End','Left Target', 'Right Target', 'Left Foot COM', 'Location', 'Best')
axis([ZMIN ZMAX YMIN YMAX])
hold off
end

function plotJointAngles(jointAngs)
% convert joint angle data from radian to deg:
%jointAngs = jointAngs*180/pi;

% Joint names:
jointName = [{'LeftHipML'}, {'LeftHipAP'}, {'RightHipML'}, {'RightHipAP'}, ...
    {'LeftAnkleML'}, {'LeftAnkleAP'}, {'RightAnkleML'}, {'RightAnkleAP'},...
    {'LeftKneeML'}, {'LeftKneeAP'}, {'RightKneeML'}, {'RightKneeAP'}];

joint_limits = [-38.8 30.5;...
                -19 121;...
                -27.75 27.75;...
                -15.3 39.7;...
                -0.1 0.1;...
                -132 0.1];
            
figTitle={'Anterior-Posterior Direction', 'MedioLateral Direction'};

% create figure with AP angles:
fID(1) = figure('units','normalized','outerposition',[0 0 1 1]);
fID(2) = figure('units','normalized','outerposition',[0 0 1 1]);
%colors = {'r','y','g','b','m','k'};

for k=1:(length(jointName)/6)
    figure(fID(k));
    title(figTitle{k});
    for j=1:(length(jointName)/4) % scan three dif joints
        subplot(1,3,j)
        plot(jointAngs(k+(j-1)*4,:),'r')%ML 1,5,9 \\ AP 2,6,10
        hold on
        plot(jointAngs(k+2+(j-1)*4,:),'b')%ML 3,7,11 \\ AP 4,8,12
        hold off
        grid on
        xlabel('Time Step')
        ylabel('Joint Angle, deg')
        xl = xlim;
        line([xl(1) xl(2)],[joint_limits(k+(j-1)*2,1) joint_limits(k+(j-1)*2,1)],'Color','r')% ML 1,3,5 \\ AP 2,4,6
        line([xl(1) xl(2)],[joint_limits(k+(j-1)*2,2) joint_limits(k+(j-1)*2,2)],'Color','r')
        legend([jointName(k+(j-1)*4), jointName(k+1+(j-1)*4)],'Location', 'Best')
    end
end
        
% for i = 1:size(jointAngs,1)
%     disp(['Using ' num2str(i) '-th row of data. Figure ID =']) %num2str(fID(mod(i+1,2)+1))])
%     disp(['fID argument is = ' num2str(mod(i+1,2)+1)])
%     figure(fID(mod(i+1,2)+1));
%     hold on
%     plot(jointAngs(i,:),colors{ceil(i/2)})
% end
% 
% for i=1:2
%    figure(fID(i));
%    grid on
%    title(figTitle{i})
%    legend([jointName(1+i-1), jointName(3+i-1), jointName(5+i-1),...
%        jointName(7+i-1), jointName(9+i-1), jointName(11+i-1)],...
%        'Location', 'Best')
%    xlabel('Time Step')
%    ylabel('Joint Angle, deg')
   %
end

        


