%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to Run Ode Solver using the blebSolver
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE

close all
clear all

% Set paramets values

R = .5;          % Initial active receptors 
alpha1 = 1/10;   % GBgamma
alpha2 = 1/12;   % GBPC
alpha3 = 1/9;    % MCOR
alpha4 = 1/16;   % RasB
alpha5 = 1/11;   % MHCKA

k1 = 1/120;    % GBgamma  
k2 = 1/130;    % GBPC
k3 = 1/145;    % MCOR
k4 = 1/160;    % RasB
k5 = 1/115;    % MHCKA

savefigure = 1; % Set 1 if want to save figure and set 0 if do not want to save figure

% Put parameter values into vector

params = [alpha1, alpha2, alpha3, alpha4, alpha5, ...
    k1, k2, k3, k4, k5];

% Set time span

tspan = [0 120];

% Set initial conditions

init_cond = [.5 0 0 0 0];

% Run Ode Solver

[T,y] = ode45(@(t,Y) blebSolver(t,Y,R,params) , tspan ,...
    init_cond);

% Steady State Solutions

VarNames = {'GBG Steady State', 'GBPC Steady State', 'MCOR Steady State',...
    'RASB Steady State', 'MHCKA Steady State'};
steady_state_table = table(y(end,1),y(end,2),y(end,3),y(end,4),y(end,5), 'VariableNames',VarNames)

% Plot solutions

figure(1)

plot(T,y, 'linewidth', 2)

title("Time Versus Concentrations", 'FontSize', 20)
xlabel("Time (Seconds)",'FontSize', 17)
ylabel("Concentrations",'FontSize', 17)

legend('GBG', 'GBPC', 'MCOR', 'RasB', 'MHCKA','location','southeast')

% Save Time Series Plot as JPG File in a Folder with the
% Date

if savefigure == 1
    
    DateDay = datestr(now, 'dd-mmm-yyyy'); % Get current date
    
    DateTime = datestr(now, 'HH:MM:SS'); % Get current time
    
    if ~exist(DateDay, 'dir')
        
       mkdir(DateDay) % Make a Directory with the Current Date if it does not already exist
       
    end
    
    fileName = strcat('/Figure', DateTime ,'.jpg'); % Name figure file name based on current time
    
    dirPath = strcat('/',DateDay, fileName); % Directory Path
    
    saveas(figure(1),[pwd dirPath]); % Save Figure in Folder
  

end