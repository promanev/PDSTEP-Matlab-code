  function cmaes  
  % (mu/mu_w, lambda)-CMA-ES 
  % CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
  % nonlinear function minimization. To be used under the terms of the
  % GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
  % Copyright: Nikolaus Hansen, 2003-09. 
  % e-mail: hansen[at]lri.fr
  %
  % This code is an excerpt from cmaes.m and implements the key parts
  % of the algorithm. It is intendend to be used for READING and
  % UNDERSTANDING the basic flow and all details of the CMA-ES
  % *algorithm*. Use the cmaes.m code to run serious simulations: it
  % is longer, but offers restarts, far better termination options, 
  % and supposedly quite useful output.
  %
  % URL: http://www.lri.fr/~hansen/purecmaes.m
  % References: See end of file. Last change: October, 21, 2010
  
  % Modified by B. Tries and R. Popov on December 2nd, 2013, Department of 
  % Computer Science and College of Medicine, University of Vermont, USA.
  % Class Project "Applying Evolutionary Strategies to Robot Locomotion" 
  % for Evolutionary Computation CS352 by Prof. Margaret J. Eppstein.
  % Corresponding author: rpopov@uvm.edu.
  
 

  % --------------------  Initialization --------------------------------  
  % User defined input parameters (need to be edited)
  
startTime = tic; %Begin timer

% PDSTEP_demo params:
num_input = 2;
num_hidden = 5;
num_output = 8;
wts_size = [num_input, num_hidden, num_output];
N = num_input*num_hidden + num_hidden*num_output;

xmean = rand(N,1);    % objective variables initial point
sigma = 0.5;          % coordinate wise standard deviation (step size)
%stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
%CHANGED: we stop only when time runs out
stopfitness = 50; 
stopeval = 1e3*N^2;   % stop after stopeval number of function 
  %evaluations
  
% other misc params:
EXP_ID = 'TARGETS 2';
if exist('EXP_ID','var')==0
    EXP_ID = 'TEST';
end
timeStamp = fix(clock);
folderName = [EXP_ID ' ' mat2str(timeStamp)];

  
% Strategy parameter setting: Selection  
lambda = 4+floor(3*log(N));  % population size, offspring number
disp(['Population size: ' num2str(lambda)])
mu = lambda/2;               % number of parents/points for recombination

weights = log(mu+1/2)-log(1:mu)';%muXone array for weighted recombination
mu = floor(mu);        
weights = weights/sum(weights);% normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2);%variance-effectiveness ofsumw_i x_i

 
% Strategy parameter setting: Adaptation
cc = (4 + mueff/N) / (N+4 + 2*mueff/N);%time constant for cumulation for 
cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));%andfor rank-mu
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma 
                                                      % usually close to 1
                                                      
% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
B = eye(N,N);                       % B defines the coordinate system
D = ones(N,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2 
eigeneval = 0;                      % track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of 
                                      %   ||N(0,I)|| == norm(randn(N,1))
%CHANGED: out.dat = []; out.datx = [];% for plotting output
max_gen = floor(stopeval/lambda);

% to cap max_gen:
if max_gen > 1000
    max_gen = 1000;
end

disp(['Starting evolution with ' num2str(max_gen) ' generations.'])
plotFitness = zeros(1,max_gen);
plotCounts = zeros(1,max_gen);
genCount = 0;

  % -------------------- Generation Loop --------------------------------
counteval = 0;  % evaluation counter to keep track of how many passed, not used currently to stop evolution
while genCount < max_gen
    % for each member in population:
    startGenTime = tic;
  
    % Generate and evaluate lambda offspring
    for k=1:lambda
      arx(:,k) = xmean + sigma * B * (D.*randn(N,1));%m+sig*Normal(0,C)
      currPopMember = arx(:,k);
      store_weights(currPopMember,wts_size);
      %arfitness(k) = feval(strfitnessfct, arx(:,k));%objectivefunctioncall
      arfitness(k) = read_fit;
      counteval = counteval+1;
    end
    
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness,2,'descend');  
    %CHANGED: now we maximize here (used to be minimization) 
    xold = xmean;
    xbest = arx(:,arindex(1)); % to save best weights
    xmean = arx(:,arindex(1:mu)) * weights;  %recombination, new mean value

    % Cumulation: Update evolution paths
    ps = (1-cs) * ps ... 
          + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
    pc = (1-cc) * pc ...
          + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma; 

    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  
    % mu difference vectors
    C = (1-c1-cmu) * C ...                   % regard old matrix  
         + c1 * (pc * pc' ...                % plus rank one update
                 + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
          + cmu * artmp * diag(weights) * artmp'; % plus rank mu update 

    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
    
    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
      eigeneval = counteval;
      C = triu(C) + triu(C,1)'; % enforce symmetry
      [B,D] = eig(C);           % eigen decomposition, 
      %B==normalized eigenvectors
      D = sqrt(diag(D));        % D contains standard deviations now
      invsqrtC = B * diag(D.^-1) * B';
    end
    
    % Break, if fitness is good enough or condition exceeds 1e14, 
    %better termination methods are advisable 
    %CHANGED: if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
    if arfitness(1) >= stopfitness
      break;
    end

    %CHANGED:  Output 
    %CHANGED: more off;  % turn pagination off in Octave
    %CHANGED: disp([num2str(counteval) ': ' num2str(arfitness(1)) ' ' ... 
    %CHANGED:       num2str(sigma*sqrt(max(diag(C)))) ' ' ...
    %CHANGED:       num2str(max(D) / min(D))]);
    %CHANGED:  with long runs, the next line becomes time consuming
    %CHANGED: out.dat = [out.dat; arfitness(1) sigma 1e5*D' ]; 
    %CHANGED: out.datx = [out.datx; xmean'];
    genCount = genCount+1; %CHANGED 
    plotFitness(genCount)= arfitness(1); %CHANGED 
    plotCounts(genCount)= genCount; %CHANGED
    
    % evolution status update every 25% of maxGen
    if genCount/max_gen==0.25|genCount/max_gen==0.5|genCount/max_gen==0.75
        disp(['Done with ' num2str(genCount/max_gen*100) '% of generations'])
    end
    
    disp(['Gen: ' int2str(genCount) ' out of ' int2str(max_gen) ', best fit: ' num2str(arfitness(1))])
    stopGenTime = toc(startGenTime);
    mins = floor(stopGenTime/60);
    secs = round(mod(stopGenTime,60));
    disp(['Generation time: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])
    
    % save best every 10 generations (useful to have a partial solution in case of early stop by user)
    if mod(genCount,10)==0
        save_best(xbest,wts_size,folderName);
    end
    
end % while, end generation loop

% plot the best fitness (taken from every generation)
plot(plotCounts, plotFitness)
% save fitness plot:
saveas(gcf,[EXP_ID ' best fitness plot'],'png');
movefile([EXP_ID ' best fitness plot.png'],folderName);

% end status 
disp(['CMAES finished on ' num2str(genCount) ' generation with fitness value of ' num2str(arfitness(1))])
stop = toc(startTime);
mins = floor(stop/60);
secs = round(mod(stop,60));
disp(['Script time: ' num2str(mins) ' minutes ' num2str(secs) ' seconds.'])

% run the best robot (DEMO version with graphics)
disp('...Starting DEMO robot')
save_best(xbest,wts_size,folderName); % saves the best pop'n member
store_weights(xbest,wts_size); % creates a weights file for running demo
system('PDSTEP_demo.exe');


%%Plotting
%Turned off for general script
%plot(plotFitness,'k-');
%xlabel('Generations');
%ylabel('Fitness');
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
    error('!!! MISSING GLUT32.DLL !!!')
end

if exist([execName '.exe'],'file')==0
    error(['!!! MISSING ' execName '.exe !!!'])
end

if exist([execName '.ilk'],'file')==0
    error(['!!! MISSING ' execName '.ilk !!!'])
end

if exist([execName '.pdb'],'file')==0
    error(['!!! MISSING ' execName '.pdb !!!'])
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
    error('!!! MISSING GLUT32.DLL !!!')
end

if exist([demoName '.exe'],'file')==0
    error(['!!! MISSING ' demoName '.exe !!!'])
end

if exist([demoName '.ilk'],'file')==0
    error(['!!! MISSING ' demoName '.ilk !!!'])
end

if exist([demoName '.pdb'],'file')==0
    error(['!!! MISSING ' demoName '.pdb !!!'])
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

% function f=robot(x)
% %Custom made fitness function. Saves current member of population in form
% %of a 4-by-8 matrix in a weights.txt file and starts the robot simulation
% %application. After that the Matlab waits for the simulation application 
% %to produce the fits.txt containing the value of the fitness of the 
% %current run. 
% 
% %%Create the weights.txt file
% fid = fopen('weights.txt','w+');
% for j=1:16
%     if mod(j,8)==0
%         fprintf(fid,'%2.7f\r\n',x(j));
%     else
%         fprintf(fid,'%2.7f ',x(j));
%     end
% end
% fclose(fid);
% %%Call the robot app
% system('final.exe');
% while (exist('APA.txt','file')==0)
%   pause(0.2);%Wait until the fits.txt appears.
% end
% %%Read the fitness value of the run
% 
% %Read the Z-traveled distance
% fid = fopen('APA.txt');
% f = fscanf(fid,'%f');
% fclose(fid);
% %disp(['Current fitness is: ' num2str(f)]);
% delete('APA.txt');
% 
% function save_best(x)
% %Function to save the best weights in the population. The logic is similar
% %to fitness function, but this function is applied to the best guy in the
% %population only. 
% fid = fopen('bestweights.txt','w+');
% for j=1:16
%     if mod(j,8)==0
%         fprintf(fid,'%2.7f\r\n',x(j));
%     else
%         fprintf(fid,'%2.7f ',x(j));
%     end
% end
% fclose(fid);

% ---------------------------------------------------------------  
%%% REFERENCES
%
% Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% Strategy on Multimodal Test Functions.  Eighth International
% Conference on Parallel Problem Solving from Nature PPSN VIII,
% Proceedings, pp. 282-291, Berlin: Springer. 
% (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
% 
% Further references:
% Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% 9(2), pp. 159-195.
% (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
%
% Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% Time Complexity of the Derandomized Evolution Strategy with
% Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
%

