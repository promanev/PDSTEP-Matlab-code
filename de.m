function de()
% d - # of dimensions in the problem to solve
% np - # of vectors in population. rule of thumb: from 5*d to 10*d, but not
% less than 4, so that proper mutation and crossover can be performed
% f - scaling factor, 0.5
% cr - crossover criterion, usually 0.1, but 0.9 and 1 can be attempted to
% see if the problem can be solved easy.

startTime = tic;

% clear prevoius graphs:
close all

% Load/init settings:
if exist('DEsettings.mat','file') ~= 0
    load('DEsettings.mat');
    disp('...Using stored DE settings: ')
    fdNames = fieldnames(DEsettings);
    for i = 1:numel(fdNames)
        disp([fdNames{i} ' = ' num2str(DEsettings.(fdNames{i}))]);
    end
    key = input('Are these settings OK? (y/n): ', 's');
    if key~='y'
        error('Initialization stopped. Change DE settings and save as DEsettings.m')
    end
else
    disp('...Using DEFAULT DE settings: ')
    DEsettings.numInput = 2;
    DEsettings.numHidden = 2;
    DEsettings.numOutput = 12;
    DEsettings.DEmethod = 0;
    DEsettings.fitFileName = 'fit.txt';
    DEsettings.demoFileName = 'PDSTEP_demo'; %.exe is added later on
    DEsettings.execFileName = 'PDSTEP_train';
    DEsettings.np = 5*(DEsettings.numInput*DEsettings.numHidden +...
        DEsettings.numHidden*DEsettings.numOutput);
    DEsettings.f = 0.5;
    DEsettings.cr = 0.1;
    DEsettings.expID = ['TWO-TARG KNEES DE2 NF-JFIX-UP-PH POP' num2str(DEsettings.np)];
    DEsettings.useScaffolding = 0;
    DEsettings.seedElite = 0;
    DEsettings.useElite = 0;
    DEsettings.maxGen = 500;
    DEsettings.maxFit = 200.0;
    DEsettings.maxStep = 400;
    DEsettings.evolveTau = 1;
    DEsettings.numTau = DEsettings.numInput + DEsettings.numHidden...
        + DEsettings.numOutput;
    DEsettings.tauRange = [0.1,20];
    DEsettings.tauMutRate = 0.10;
    
    % NF - new fitness, JFIX - Josh's fix for forward reaching
    % UP - simulations where robots that fall below a threshold are
    % terminated, PH - pelvis height term is calculated by the same formula
    % as forward reaching term. 
    
    % Make sure that these settings are correct:
    fdNames = fieldnames(DEsettings);
    for i = 1:numel(fdNames)
        disp([fdNames{i} ' = ' num2str(DEsettings.(fdNames{i}))]);
    end
    key = input('Are these settings OK? (y/n): ', 's');
    if key~='y'
        error('Initialization stopped. Change DE settings and save as DEsettings.mat')
    end
end
  
% PDSTEP_demo params:
num_input = DEsettings.numInput;
num_hidden = DEsettings.numHidden;
num_output = DEsettings.numOutput;
% weight matrices dimensions:
wts_size = [num_input, num_hidden, num_output];
% vector size:
d = num_input*num_hidden + num_hidden*num_output;

% name of file where fitness is recorded by C++ code:
fileName = DEsettings.fitFileName;
% DE configuration is either rand/1/bin (1) or best/2/bin(0):
DEmethod = DEsettings.DEmethod;
% init default parameters:
np = DEsettings.np;
f = DEsettings.f;
cr = DEsettings.cr;
max_gen = DEsettings.maxGen; % limits the generations
max_fit = DEsettings.maxFit; % limits the fitness
% check that np > 4
if np < 4
    error('Not enough population members NP; needs to be 4 at least!')
end
% Experiment ID to be used for folder naming: 
EXP_ID = DEsettings.expID;
EXP_ID = [EXP_ID ' ' num2str(max_gen) 'GEN'];
EXP_ID = [EXP_ID ' SS' num2str(DEsettings.maxStep)];
% Add a time stamp:
timeStamp = fix(clock);
folderName = [EXP_ID ' ' mat2str(timeStamp)];

% bool value for scaffolding option:
useScaffolding = DEsettings.useScaffolding;
% bool value for seeding from an elite:
seedElite = DEsettings.seedElite;
% bool value for using elite during evolution:
useElite = DEsettings.useElite;
% Add these tags to experiment ID, if these features are used:
if useScaffolding
    EXP_ID = [EXP_ID ' S'];
end
if seedElite
    EXP_ID = [EXP_ID ' sE'];
end
if useElite
    EXP_ID = [EXP_ID ' E'];
end
if DEsettings.evolveTau
    EXP_ID = [EXP_ID ' ETau'];
end
% number of elite members:
if useElite
    eliteNum = ceil(np*0.05);
    eliteArr = zeros(eliteNum,d);
    eliteFit = zeros(1,eliteNum);
end

% NOT USED:
% To use if there is necessity to terminate evolution early after
% "breakCounter" generations where fitness improvement was less than
% "min_fit_change". Not used atm, because evolution can stick to the same
% fitness for very long time. 
% min_fit_change = 0.05; % minimal change in fitness between generations necessary for continuation of evolution
% breakCounter = 0; %counts how many consecutive generations were without improvement. 

count = 1;
trial = zeros(1,d); % trial vector to be tested against target vector
if DEsettings.evolveTau
    trialTau = zeros(1,DEsettings.numTau);
    tempTau = zeros(np,DEsettings.numTau); %array for swapping
end
fitness = zeros(1,np); % vector for current fitness values
bestfit = zeros(1,max_gen); % vector for best fitness values

% seed the population:
if seedElite
    elite = read_seed;
    x = seed_pop(elite);
else
    if exist('seed.mat','file')~=0
        load seed
        rng(seed);
    end
    x = rand(np,d); % each row is a member of population, each column - dimension
    seed = rng;
end

% initialize tau array:
if DEsettings.evolveTau
    tau = DEsettings.tauRange(2)*rand(np,DEsettings.numTau);
    tau(tau<DEsettings.tauRange(1)) = DEsettings.tauRange(1);
end
    
tempx = zeros(np,d); % array for swapping

% calculate their fitness the first time:
initTime = tic;
for l=1:np
    currPopMember = x(l,:);
    store_weights(currPopMember,wts_size);
    if DEsettings.evolveTau
        currTau = tau(l,:);
        storeTau(currTau,wts_size);
    end
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
        if DEmethod
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
                
                j = mod((j),d)+1;
                
                if DEsettings.evolveTau
                    m=ceil(rand*DEsettings.numTau);
                    for n=1:DEsettings.numTau
                        if or(rand < DEsettings.tauMutRate, k==DEsettings.numTau)
                        trialTau(m) = tau(members(3),m) + f*(tau(members(1),m) - tau(members(2),m));
                        else
                        trialTau(m) = tau(i,m);
                        end

                        % re-normalize tau:
                        if trialTau(m) < DEsettings.tauRange(1)
                            trialTau(m) = DEsettings.tauRange(1);
                        end

                        if trialTau(m) > DEsettings.tauRange(2)
                            trialTau(m) = DEsettings.tauRange(2);
                        end
                        m = mod((m),DEsettings.numTau)+1; % loops j values around starting from any point in vector of all values
                    end
                end % end tau part 
            end
        else
            % MUTATE|RECOMBINE
            
            % Find the best vector in the current population:
            [~,maxInx] = max(fitness);
            currBestX = x(maxInx,:);
            if DEsettings.evolveTau
                currBestTau = tau(maxInx,:);
            end
            % randomly pick FOUR members:
            members = randperm(np,4);
            % randomly pick the dimension from which to start mutation:
            j = ceil(rand*d);
            for k=1:d

                if or(rand < cr, k==d)
                    trial(j) = currBestX(j) + f*(x(members(1),j) + x(members(2),j) - x(members(3),j) - x(members(4),j));
                else
                    trial(j) = x(i,j);
                end 
                j = mod((j),d)+1;
                % if also evolve tau:
                if DEsettings.evolveTau
                    m=ceil(rand*DEsettings.numTau);
                    for n=1:DEsettings.numTau
                        if or(rand < DEsettings.tauMutRate, k==DEsettings.numTau)
                        trialTau(m) = currBestTau(m) + f*(tau(members(1),m) + tau(members(2),m) - tau(members(3),m) - tau(members(4),m));
                        else
                        trialTau(m) = tau(i,m);
                        end
                    
                        % re-normalize tau:
                        if trialTau(m) < DEsettings.tauRange(1)
                            trialTau(m) = DEsettings.tauRange(1);
                        end
                    
                        if trialTau(m) > DEsettings.tauRange(2)
                            trialTau(m) = DEsettings.tauRange(2);
                        end
                        m = mod((m),DEsettings.numTau)+1; % loops j values around starting from any point in vector of all values
                    end
                end % end tau evolve
                
            end
        end
        
%         % Mutate taus:
%         if DEsettings.evolveTau
%             tempTau = tau;
%            for l=1:length(tau)
%               if rand<DEsettings.tauMutRate
%                   tempTau(l) = DEsettings.tauRange(2)*rand;
%                   if tempTau(l)<DEsettings.tauRange(1)
%                       tempTau(l) = DEsettings.tauRange(1);
%                   end
%               end
%            end
%         end
        % EVALUATE|SELECT:
        if DEsettings.evolveTau
            storeTau(trialTau, wts_size);
        end
        store_weights(trial,wts_size);
        
        % choose appropriate training robot file
        if useScaffolding
            if (count<max_gen*0.75)
                execName = DEsettings.execFileName; % no upright punishment
            else
                execName = 'PDSTEP_train_up';
            end
        else
            execName = DEsettings.execFileName;
        end
        
        score = read_fit(fileName,execName);
        
        if score>=fitness(i)%maximize fitness!
            tempx(i,:) = trial(:);
            if DEsettings.evolveTau
                tempTau(i,:) = trialTau(:);
            end
            fitness(i) = score;
        else
            tempx(i,:) = x(i,:);
            if DEsettings.evolveTau
                tempTau(i,:) = tau(i,:);
            end
        end
        
    end % end population loop
    
    % SWAP temp array and population array:
    x = tempx;
    if DEsettings.evolveTau
        tau = tempTau;
    end
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
    disp(['Gen: ' int2str(count) ' out of ' int2str(max_gen) ', best fit: ' num2str(bestfit(count))])
    count = count + 1;
    
    stopGenTime = toc(startGenTime);
    mins = floor(stopGenTime/60);
    secs = round(stopGenTime - mins*60);
    if secs < 0
        secs = 0;
    end
    disp(['Generation time: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])
    
    timeLeft = stopGenTime*max_gen - count*stopGenTime;
    hours = floor(timeLeft/3600);
    if hours < 0
        hours = 0;
    end
    mins = floor((timeLeft - hours*3600)/60);
    if mins < 0
        mins = 0;
    end
    secs = round(timeLeft - hours*3600 - mins*60);
    if secs < 0
        secs = 0;
    end
    disp(['Approx. time left: ' num2str(hours) ' hours ' num2str(mins) ' minutes ' num2str(secs) ' seconds until all ' num2str(max_gen) ' generations are done.'])
    
    % save best every 10 generations (useful to have a partial solution in case of early stop by user)
    if mod(count,10)==0
        [~,inx] = max(fitness);
        save_best(x(inx,:),wts_size,folderName,['weights' num2str(count) 'gen']);
        if DEsettings.evolveTau
            saveBestTau(tau(inx,:),wts_size,folderName,['tau' num2str(count) 'gen']);
        end
        save([folderName '\' EXP_ID]);% Saves all of the variable from the run in the EXP_ID.mat
        % Do this only once
        if count == 10
            % Saving seed/settings used for this run:
            root = cd(folderName);
            save('seed','seed');
            save('DEsettings', 'DEsettings');
            cd(root);
            % keep a copy of settings to be used for next runs:
            save('DEsettings', 'DEsettings');
            
            % move source files to the folder:
            if exist([folderName '\source'],'dir')==0
                disp(['Folder ' folderName '\source does not exist. Creating this folder...'])
                mkdir([folderName '\source']);
            end
            oldLocation = cd([folderName '\source']);
            if exist('PDSTEP_demo.cpp','file')==0
                cd(oldLocation);
                disp('Moving source files...')
                copyfile('C:\Users\Roman\Documents\Visual Studio 2015\Projects\bullet-2.82-r2704\Demos\PDSTEP_demo\PDSTEP_demo.cpp',[folderName '\source'],'f');
                copyfile('C:\Users\Roman\Documents\Visual Studio 2015\Projects\bullet-2.82-r2704\Demos\PDSTEP_demo\PDSTEP_demo.h',[folderName '\source'],'f');
                copyfile('C:\Users\Roman\Documents\Visual Studio 2015\Projects\bullet-2.82-r2704\Demos\PDSTEP_demo\main.cpp',[folderName '\source'],'f');
            else
                cd(oldLocation);
            end %end move source files 
        end %end saving procedure done on 10-th turn
    end %end backing up done every 10 turns 
end %end generational loop

% plot the best fitness (taken from every generation)
if length(1:count-1)==length(bestfit(bestfit>0))
    plot(1:count-1, bestfit(bestfit>0))
else
    warning('!!!count vector and bestfit are not equal, plotting only bestfit!!!')
    plot(bestfit(bestfit>0))
end
xlabel('Generations'); ylabel('Fitness')
title(['Results for "' EXP_ID '" experiment'])
% save fitness plot:
saveas(gcf,[EXP_ID ' best fitness plot'],'png');
movefile([EXP_ID ' best fitness plot.png'],folderName);

% move any readme files to experiment's folder:
readmeList = dir('README*.txt');
if size(readmeList,1)==1
    movefile(readmeList.name,folderName);
elseif size(readmeList,1)>1
    for k=1:size(readmeList,1)
         movefile(readmeList(k).name,folderName);
    end
end

% end status 
disp(['DE finished on ' num2str(count-1) ' generation with fitness value of ' num2str(max(fitness))])
stop = toc(startTime);
hours = floor(stop/3600);
if hours < 0
    hours = 0;
end
mins = floor((stop - hours*3600)/60);
if mins < 0
    mins = 0;
end
secs = round(stop - hours*3600 - mins*60);
if secs < 0
    secs = 0;
end
disp(['Script time: ' num2str(hours) ' hours ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])

% run the best robot (DEMO version with graphics)
disp('...Starting DEMO robot')
[~,inx] = max(fitness);
save_best(x(inx,:),wts_size,folderName); % saves the best pop'n member
if DEsettings.evolveTau
    saveBestTau(tau(inx,:),wts_size,folderName); % saves the best tau member
end
store_weights(x(inx,:),wts_size); % creates a weights file for running demo
cd(folderName);
system('PDSTEP_demo.exe');
end

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

end

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
end

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
end

function save_best(wts_vector,wts_size,folderName, fileName)
% create a default folder, if no folder name is supplied
if nargin<4
    fileName = 'weights';
end

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
oldFolder = cd(folderName);
if ~exist('PDSTEP*','file')
    cd(oldFolder);
    copyfile('PDSTEP_demo.exe', folderName);
    copyfile('PDSTEP_demo.ilk', folderName);
    copyfile('PDSTEP_demo.pdb', folderName);
    copyfile('PDSTEP_train.exe', folderName);
    copyfile('PDSTEP_train.ilk', folderName);
    copyfile('PDSTEP_train.pdb', folderName);
    copyfile('GLUT32.DLL', folderName);
else
    cd(oldFolder);
end

oldFolder = cd(folderName);

fid = fopen([fileName '.txt'],'w+');
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
end

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
end

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
end

function storeTau(tau,wts_size)
% Store tau array in a text file:
fid = fopen('tau.txt','w+');

for j=1:length(tau) 
    if j<wts_size(1)
        fprintf(fid,'%2.4f ',tau(j));
    elseif j==wts_size(1)
        fprintf(fid,'%2.4f\r\n',tau(j));
    elseif and(j>wts_size(1), j<(wts_size(1) + wts_size(2)))
        fprintf(fid,'%2.4f ',tau(j));
    elseif j==(wts_size(1) + wts_size(2))
        fprintf(fid,'%2.4f\r\n',tau(j));
    elseif j>(wts_size(1) + wts_size(2))
        fprintf(fid,'%2.4f ',tau(j));
    end
end
fclose(fid);
end

function saveBestTau(tau,wts_size,folderName, fileName)
if nargin < 4
    fileName = 'tau';
end
fileName = [fileName '.txt'];
% Store tau array in a text file:
oldFolder = cd(folderName);
fid = fopen(fileName,'w+');
for j=1:length(tau) 
    if j<wts_size(1)
        fprintf(fid,'%2.4f ',tau(j));
    elseif j==wts_size(1)
        fprintf(fid,'%2.4f\r\n',tau(j));
    elseif and(j>wts_size(1), j<(wts_size(1) + wts_size(2)))
        fprintf(fid,'%2.4f ',tau(j));
    elseif j==(wts_size(1) + wts_size(2))
        fprintf(fid,'%2.4f\r\n',tau(j));
    elseif j>(wts_size(1) + wts_size(2))
        fprintf(fid,'%2.4f ',tau(j));
    end
end
cd(oldFolder);
fclose(fid);
end