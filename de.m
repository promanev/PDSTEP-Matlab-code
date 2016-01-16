function de(d,np,f,cr)
% d - # of dimensions in the problem to solve
% np - # of vectors in population. rule of thumb: from 5*d to 10*d, but not
% less than 4, so that proper mutation and crossover can be performed
% f - scaling factor, 0.5
% cr - crossover criterion, usually 0.1, but 0.9 and 1 can be attempted to
% see if the problem can be solved easy.

startTime = tic;

% clear prevoius graphs:
close all

% PDSTEP_demo params:
num_input = 2;
num_hidden = 2;
num_output = 8;
wts_size = [num_input, num_hidden, num_output];

if nargin<1
    d = num_input*num_hidden + num_hidden*num_output;
end
% name of file where fitness is recorded by C++ code:
fileName = 'fit.txt';

% init default parameters:
if nargin<2
    np = 50;
    f = 0.5;
    cr = 0.1;
    %disp(['Using default DE params: 1) NP = ' num2str(np) ', F = ' num2str(f) ', CR = ' num2str(cr)])
end

EXP_ID = ['DEBUG DE POP' num2str(np)];

% check that np > 4
if np < 4
    error('Not enough population members NP; needs to be 4 at least!')
end

% bool value for scaffolding option:
useScaffolding = 0;
if useScaffolding
    EXP_ID = [EXP_ID ' S'];
else
    EXP_ID = [EXP_ID ' NS'];
end

% bool value for seeding from an elite:
seedElite = 0;
if seedElite
    EXP_ID = [EXP_ID ' sE'];
else
    EXP_ID = [EXP_ID ' NsE'];
end

% bool value for using elite during evolution:
useElite = 0;
if useElite
    EXP_ID = [EXP_ID ' E'];
else
    EXP_ID = [EXP_ID ' NE'];
end

% number of elite members:
if useElite
    eliteNum = ceil(np*0.05);
    eliteArr = zeros(eliteNum,d);
    eliteFit = zeros(1,eliteNum);
end

max_gen = 200; % limits the generations
max_fit = 10^3; % limits the fitness
EXP_ID = [EXP_ID ' ' num2str(max_gen) 'GEN'];

% naming scheme (add. stuff, not yet automatized):
EXP_ID = [EXP_ID ' SS100 TXT'];


if exist('EXP_ID','var')==0
    EXP_ID = 'TEST';
end
timeStamp = fix(clock);
folderName = [EXP_ID ' ' mat2str(timeStamp)];



% To use if there is necessity to terminate evolution early after
% "breakCounter" generations where fitness improvement was less than
% "min_fit_change". Not used atm, because evolution can stick to the same
% fitness for very long time. 
% min_fit_change = 0.05; % minimal change in fitness between generations necessary for continuation of evolution
% breakCounter = 0; %counts how many consecutive generations were without improvement. 

count = 1;
trial = zeros(1,d); % trial vector to be tested against target vector
fitness = zeros(1,np); % vector for fitness values
bestfit = zeros(1,max_gen);

% seed the population:
if seedElite
    elite = read_seed;
    x = seed_pop(elite);
else
    x = rand(np,d); % each row is a member of population, each column - dimension
end

tempx = zeros(np,d); % array for swapping

% calculate their fitness the first time:
initTime = tic;
for l=1:np
    currPopMember = x(l,:);
    store_weights(currPopMember,wts_size);
    fitness(l) = read_fit;
end
stopInitTime = toc(initTime);
mins = floor(stopInitTime/60);
secs = round(mod(stopInitTime,60));
disp(['Initialization of pop.array took: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])
disp(['Or ' num2str(stopInitTime/np) ' seconds per pop.member'])

% if using elite, populate the elite matrix:
if useElite
    temp = fitness;
    for m=1:eliteNum
       [~,index] = max(temp);
       eliteArr(m,:) = x(index,:);
       eliteFit(m) = temp(index);
       % exclude this fitness value from future consideration:
       temp(index) = 0;
    end
    % seeding from elite, use it as the last member of elite array:
    if seedElite
        eliteArr(eliteNum,:) = read_seed;
        copyfile('seed.txt','weights.txt');
        eliteFit(eliteNum) = read_fit;
    end
    clear temp
    clear index
end


% Main loop:
while and(count <= max_gen, min(fitness)<max_fit)
    % for each member in population:
    startGenTime = tic;
    % check if previous fitness increase was smaller than threshold
%     if count>1
%         if breakCounter > 10
%             disp(['Evolution stopped due to lack of fitness improvement at ' int2str(count-1) ' generation.'])
%             break
%         end
%         
%         if (bestfit(count)-bestfit(count-1))<min_fit_change
%             breakCounter = breakCounter + 1;
%         end
%     end

    % go through each population member:
    for i = 1:np
        % MUTATE|RECOMBINE
        % randomly pick three members:
        members = randperm(np,3);
        % randomly pick the dimension from which to start mutation:
        j = ceil(rand*d);
        for k=1:d
            
            if or(rand < cr, k==d)
                trial(j) = x(members(3),j) + f*(x(members(1),j) - x(members(2),j));
            else
                trial(j) = x(i,j);
            end        
            j = mod((j),d)+1; % loops j values around starting from any point in vector of all values
        end
        
        % EVALUATE|SELECT:        
        store_weights(trial,wts_size);
        
        % choose appropriate training robot file
        if useScaffolding
            if (count<max_gen*0.7)
                execName = 'PDSTEP_train_nu'; % no upright punishment
            else
                execName = 'PDSTEP_train';
            end
        else
            execName = 'PDSTEP_train';
        end
        
        score = read_fit(fileName,execName);
        
        if score>=fitness(i)%maximize fitness!
            tempx(i,:) = trial(:);
            fitness(i) = score;
        else
            tempx(i,:) = x(i,:);
        end
        
    end % end population loop
    
    % SWAP temp array and population array:
    x = tempx;
    
    % update elite array and/or replace members of population with elites:
    if useElite
        temp = fitness;
        % for each elite, check if it needs to be replaced
        for n=1:eliteNum
            [~,index] = max(temp);
            if eliteFit(n)<temp(index)
                eliteArr(n,:) = x(index,:);
                eliteFit(n) = fitness(index);
            % if not, then replace the weakest member of population with this elite    
            else
                [~,minIndex] = min(fitness);
                x(minIndex,:) = eliteArr(n,:);
                fitness(minIndex) = eliteFit(n);
            end
            % to exclude this population member, equal its temp. fitness
            % value to a zero. 
            temp(index) = 0;
        end
        
    end
    % evolution status update every 25% of maxGen
    if count/max_gen==0.25|count/max_gen==0.5|count/max_gen==0.75
        disp(['Done with ' num2str(count/max_gen*100) '% of generations'])
    end
    
    bestfit(count) = max(fitness);
%     disp(['Best fitness this gen: ' num2str(max(fitness)) ', and it is recorded in bestfit(' num2str(count) ') = ' num2str(bestfit(count))])
    disp(['Gen: ' int2str(count) ' out of ' int2str(max_gen) ', best fit: ' num2str(bestfit(count))])
    count = count + 1;
    stopGenTime = toc(startGenTime);
    mins = floor(stopGenTime/60);
    secs = round(mod(stopGenTime,60));
    disp(['Generation time: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])
    timeLeft = stopGenTime*max_gen - count*stopGenTime;
    hours = floor(timeLeft/3600);
    mins = floor((timeLeft - hours*3600)/60);
    secs = round(timeLeft - hours*3600 - mins*60);
    disp(['Approx. time left: ' num2str(hours) ' hours ' num2str(mins) ' minutes ' num2str(secs) ' seconds until all ' num2str(max_gen) ' generations are done.'])
    
    % save best every 10 generations (useful to have a partial solution in case of early stop by user)
    if mod(count,10)==0
        [~,inx] = max(fitness);
        save_best(x(inx,:),wts_size,folderName);
        % move source files to the folder:
        if exist([folderName '\source'],'dir')==0
            mkdir([folderName '\source']);
        end
        oldLocation = cd([folderName '\source']);
        if exist('PDSTEP_demo.cpp','file')==0
            cd(oldLocation);
            copyfile('C:\Users\Roman\Documents\Visual Studio 2015\Projects\bullet-2.82-r2704\Demos\PDSTEP_demo\PDSTEP_demo.cpp',[folderName '\source'],'f');
            copyfile('C:\Users\Roman\Documents\Visual Studio 2015\Projects\bullet-2.82-r2704\Demos\PDSTEP_demo\PDSTEP_demo.h',[folderName '\source'],'f');
            copyfile('C:\Users\Roman\Documents\Visual Studio 2015\Projects\bullet-2.82-r2704\Demos\PDSTEP_demo\main.cpp',[folderName '\source'],'f');
        else
            cd(oldLocation);
        end
        
    end
end
% plot the best fitness (taken from every generation)
if length(1:count)==length(bestfit(bestfit>0))
    plot(1:count, bestfit(bestfit>0))
else
    warning('!!!count vector and bestfit are not equal, plotting only bestfit!!!')
    plot(bestfit(bestfit>0))
end
xlabel('Generations'); ylabel('Fitness (max = 1)')
title(['Evolution results for "' EXP_ID '" experiment'])
% save fitness plot:
saveas(gcf,[EXP_ID ' best fitness plot'],'png');
movefile([EXP_ID ' best fitness plot.png'],folderName);

% end status 
disp(['DE finished on ' num2str(count) ' generation with fitness value of ' num2str(max(fitness))])
stop = toc(startTime);
hours = floor(stop/3600);
mins = floor((stop - hours*3600)/60);
secs = round(stop - hours*3600 - mins*60);
disp(['Script time: ' num2str(hours) ' hours ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])

% run the best robot (DEMO version with graphics)
disp('...Starting DEMO robot')
[~,inx] = max(fitness);
% DEBUG:
disp('Fitness vector: ')
disp(fitness)
disp(['Max fitness is at ' num2str(inx) ' vector'])
disp('Checking that the best is written to output dir')
disp(x(inx,:))
save_best(x(inx,:),wts_size,folderName); % saves the best pop'n member
store_weights(x(inx,:),wts_size); % creates a weights file for running demo
system('PDSTEP_demo.exe');

%% Misc functions:
function store_weights(wts_vector,wts_size)
% Function that stores wts in the local directory for C++ .exe to use
% during the simulation:

fid = fopen('weights.txt','w+');
wIJlength = wts_size(1)*wts_size(2);
wJKlength = wts_size(2)*wts_size(3);

for j=1:wIJlength 
    if mod(j,wts_size(2))==0
        fprintf(fid,'%2.4f\r\n',wts_vector(j));
    else
        fprintf(fid,'%2.4f ',wts_vector(j));
    end
end
    
for j=(wIJlength+1):(wJKlength + wIJlength)
    if mod((j - wIJlength),wts_size(3))==0
        fprintf(fid,'%2.4f\r\n',wts_vector(j));
    else
        fprintf(fid,'%2.4f ',wts_vector(j));
    end   
end

fclose(fid);

function f = read_fit(fileName,execName)
% Function runs the simulation and reads fitness from a FILE:
if nargin<2
    execName = 'PDSTEP_train';
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

if nargin<1
    fileName = 'fit.txt';
end

% Making sure there are no un-closed files:
if ~isempty(fopen('all'))
    fclose('all');
end

fid = fopen(fileName,'r');
f = fscanf(fid,'%f',[1,1]);
fclose(fid);
% disp(['Current fitness is: ' num2str(f)]);
delete(fileName);

function message = socket_fit(execName, numRetries)
% Function that runs simulation and receives fitness via socket:
if nargin<2
    numRetries = 30;
end

if nargin<1
    execName = 'PDSTEP_train';
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

% using JAVA IO libraries 
import java.net.ServerSocket
import java.io.*
% port is hardcoded in both C++ and here
output_port = 27015;
% initializing:
server_socket  = [];
input_socket  = [];
retry = 0;

while true
    
    retry = retry + 1;
    if ((numRetries > 0) && (retry > numRetries))
        error('Couldn''t obtain fitness!');
        break;
    end
    
    try
        server_socket = ServerSocket(output_port);
        server_socket.setSoTimeout(500);
        % run the simulation app in the background, so that Matlab continues
        % working (crucial for receiving data through the socket!)
        dos([execName '.exe']);
%         system([execName '.exe &']);
        input_socket = server_socket.accept;
        input_stream = input_socket.getInputStream;
        d_input_stream = DataInputStream(input_stream);
        message = zeros(1, 8, 'uint8');
%         disp(['Message length(before conversion): ' int2str(length(message))])
        for i = 1:8
            message(i) = d_input_stream.read;
        end
%         disp(['Got this value from socket', message])
        message = char(message);
%         disp(['Same value in char format: ', message])
        message = str2double(message);
%         disp(['Message length(after conversion): ' int2str(length(message))])
%         disp(['Same value in double format: ', num2str(message)])
        server_socket.close;
        input_socket.close;
        
        if length(message)==1
            break;
        end
        
    catch
        if ~isempty(server_socket)
            server_socket.close
        end

        if ~isempty(input_socket)
            input_socket.close
        end
        % pause before retrying
        pause(0.1);
    end
end

function save_best(wts_vector,wts_size,folderName)
% create a default folder, if no folder name is supplied
if nargin<3
    EXP_ID = 'TEST';
    % create a date-time stamp for folder name:
    timeStamp = fix(clock);
    folderName = [EXP_ID ' ' mat2str(timeStamp)];
    mkdir(folderName);
end
if exist(folderName,'dir')==0
    mkdir(folderName);
end

% check that all files are present:
demoName = 'PDSTEP_demo';
if exist('GLUT32.DLL','file')==0
    warning('!!! MISSING GLUT32.DLL !!!')
    input('Copy necessary file and press any key','s');
end

if exist([demoName '.exe'],'file')==0
    warning(['!!! MISSING ' demoName '.exe !!!'])
    input('Copy necessary file and press any key','s');
end

if exist([demoName '.ilk'],'file')==0
    warning(['!!! MISSING ' demoName '.ilk !!!'])
    input('Copy necessary file and press any key','s');
end

if exist([demoName '.pdb'],'file')==0
    warning(['!!! MISSING ' demoName '.pdb !!!'])
    input('Copy necessary file and press any key','s');
end

% copy demo robot files into the new folder:
copyfile('PDSTEP_demo.exe', folderName);
copyfile('PDSTEP_demo.ilk', folderName);
copyfile('PDSTEP_demo.pdb', folderName);
copyfile('PDSTEP_train.exe', folderName);
copyfile('PDSTEP_train.ilk', folderName);
copyfile('PDSTEP_train.pdb', folderName);
copyfile('GLUT32.DLL', folderName);

oldFolder = cd(folderName);

fid = fopen('weights.txt','w+');
wIJlength = wts_size(1)*wts_size(2);
wJKlength = wts_size(2)*wts_size(3);

for j=1:wIJlength 
    if mod(j,wts_size(2))==0
        fprintf(fid,'%2.4f\r\n',wts_vector(j));
    else
        fprintf(fid,'%2.4f ',wts_vector(j));
    end
end
    
for j=(wIJlength+1):(wJKlength + wIJlength)
    if mod((j - wIJlength),wts_size(3))==0
        fprintf(fid,'%2.4f\r\n',wts_vector(j));
    else
        fprintf(fid,'%2.4f ',wts_vector(j));
    end   
end
fclose(fid);
cd(oldFolder);

function a = read_seed(fileName)
% Function that uploads a seeding genome:
if nargin<1
    fileName = 'seed.txt';
end

fid = fopen(fileName);
a = [];
tline = fgetl(fid);
while ischar(tline)
    a = [a str2num(tline)];
    tline = fgetl(fid);
end

fclose(fid);

function np = seed_pop(elite, popNum, wts_size)
% function that seeds initial population from a single (need several?)
% elite member from previous evolutionary runs. 

if nargin<3
% PDSTEP_demo default params:
    num_input = 2;
    num_hidden = 2;
    num_output = 8;
    wts_size = num_input*num_hidden + num_hidden*num_output;
end

if nargin<2
    popNum = 50; % originally suggested 5x - 10x wts_size, but lower due to hardware limitation (too slow!)
end

% Rate at which original elite member is mutated:
mutRate = 0.05;

np = zeros(popNum, wts_size);

for i=1:popNum
    for j=1:wts_size;
        if (rand < mutRate)
            np(i,j) = elite(j)*rand;
        else
            np(i,j) = elite(j);
        end        
    end 
end

