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
num_hidden = 5;
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

% check that np > 4
if np < 4
    error('Not enough population members NP; needs to be 4 at least!')
end

% bool value for scaffolding option:
useScaffolding = 0;

% bool value for seeding from an elite:
seedElite = 1;

% bool value for using elite during evolution:
useElite = 1;

% number of elite members:
if useElite
    eliteNum = ceil(np*0.05);
    eliteArr = zeros(eliteNum,d);
    eliteFit = zeros(1,eliteNum);
end

% other misc params:
EXP_ID = 'TARGETS DE ELITISM';
if exist('EXP_ID','var')==0
    EXP_ID = 'TEST';
end
timeStamp = fix(clock);
folderName = [EXP_ID ' ' mat2str(timeStamp)];

max_gen = 1000; % limits the generations
max_fit = 10^3; % limits the fitness

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
for l=1:np
    currPopMember = x(l,:);
    store_weights(currPopMember,wts_size);
    fitness(l) = read_fit;
end

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
while and(count < max_gen, min(fitness)<max_fit)
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
    disp(['Gen: ' int2str(count) ' out of ' int2str(max_gen) ', best fit: ' num2str(bestfit(count))])
    count = count + 1;
    stopGenTime = toc(startGenTime);
    mins = floor(stopGenTime/60);
    secs = round(mod(stopGenTime,60));
    disp(['Generation time: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])
    
    % save best every 10 generations (useful to have a partial solution in case of early stop by user)
    if mod(count,10)==0
        [~,inx] = max(fitness);
        save_best(x(inx,:),wts_size,folderName);
    end
end
% plot the best fitness (taken from every generation)
if length(1:count-1)==length(bestfit(bestfit>0))
    plot(1:count-1, bestfit(bestfit>0))
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
disp(['DE finished on ' num2str(count) ' generation with fitness value of ' num2str(min(fitness))])
stop = toc(startTime);
mins = floor(stop/60);
secs = round(mod(stop,60));
disp(['Script time: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])

% run the best robot (DEMO version with graphics)
disp('...Starting DEMO robot')
[~,inx] = max(fitness);
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
    
for j=1:wJKlength
    if mod(j,wts_size(3))==0
        fprintf(fid,'%2.4f\r\n',wts_vector(j));
    else
        fprintf(fid,'%2.4f ',wts_vector(j));
    end   
end

fclose(fid);

function f = read_fit(fileName,execName)
% Function that stores wts in the local directory for C++ .exe to use
% during the simulation:
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
% necessary?
% while (exist(fileName,'file')==0)
%   pause(0.2);
% end

% Making sure there are no un-closed files:
if ~isempty(fopen('all'))
    fclose('all');
end

fid = fopen(fileName,'r');
f = fscanf(fid,'%f',[1,1]);
fclose(fid);
%disp(['Current fitness is: ' num2str(f)]);
delete(fileName);

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
    
for j=1:wJKlength
    if mod(j,wts_size(3))==0
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
    num_hidden = 5;
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

