% play a demo simulation:
system('PDSTEP_demo.exe');
% read COM path:
fid = fopen('com.txt','r');
f = fscanf(fid,'%f',[3,Inf]);
fclose(fid);
% mirror x-axis values, so the graph is the same as looking at robot's back
% as in simulation: 
f(1,:) = -1*f(1,:);
% read target locations:
fid = fopen('targets.txt','r');
t = fscanf(fid,'%f',[3,Inf]);
fclose(fid);
t(1,:) = -1*t(1,:);

XMIN = -3.5;
XMAX = 3.5;
YMIN = -0.5;
YMAX = 3.5;
ZMIN = -0.5;
ZMAX = 3;
length = 0.812908;
width = 0.541939;

figure
subplot(1,3,1)

plot(f(1,:),f(2,:))
hold on
plot(f(1,1),f(2,1),'ro')
plot(f(1,end),f(2,end),'ko')
plot(t(1,1), t(2,1), 'gx')
plot(t(1,2), t(2,2), 'bx')
title('COM path in XY plane (BACK)')
xlabel('X-axis (LEFT \Rightarrow RIGHT)')
ylabel('Y-axis (DOWN \Rightarrow UP)')
grid on
legend('COM path', 'Start', 'End','Left Target', 'Right Target', 'Location', 'Best')
axis([XMIN XMAX YMIN YMAX])
hold off

subplot(1,3,2)

plot(f(1,:),f(3,:))
hold on
plot(f(1,1),f(3,1),'ro')
plot(f(1,end),f(3,end),'ko')
plot(t(1,1), t(3,1), 'gx')
plot(t(1,2), t(3,2), 'bx')
title('COM path in XZ plane (BIRD''S-EYE VIEW)')
xlabel('X-axis (LEFT \Rightarrow RIGHT)')
ylabel('Z-axis (BACKWARD \Rightarrow FORWARD)')
% Draw base of support:
%back line
line([(-0.3609-width/2) (0.3609+width/2)],[-length/2 -length/2])
% front line
line([(-0.3609-width/2) (0.3609+width/2)],[length/2 length/2])
% left line
line([(-0.3609-width/2) (-0.3609-width/2)],[-length/2 length/2])
%right line
line([(0.3609+width/2) (0.3609+width/2)],[-length/2 length/2])
grid on
legend('COM path', 'Start', 'End','Left Target', 'Right Target', 'Location', 'Best')
axis([XMIN XMAX ZMIN ZMAX])
hold off

subplot(1,3,3)

plot(f(3,:),f(2,:))
hold on
plot(f(3,1),f(2,1),'ro')
plot(f(3,end),f(2,end),'ko')
plot(t(3,1), t(2,1), 'gx')
plot(t(3,2), t(2,2), 'bx')
title('COM path in ZY plane (RIGHT SIDE)')
xlabel('Z-axis (BACKWARD \Rightarrow FORWARD)')
ylabel('Y-axis (DOWN \Rightarrow UP)')
grid on
legend('COM path', 'Start', 'End','Left Target', 'Right Target', 'Location', 'Best')
axis([ZMIN ZMAX YMIN YMAX])
hold off

delete('com.txt')
delete('targets.txt')