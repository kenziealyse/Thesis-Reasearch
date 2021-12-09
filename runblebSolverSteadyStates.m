%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to Run Ode Solver for the Steady 
% state ODE and save a times series plot as a figure
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE

close all
clear all

% Set paramets values

R = 0.8;          % Initial active receptors 
k1plus = 1/10;   % GBgamma
k2plus = 1/12;   % GBPC
k3plus = 1/9;    % MCOR
k4plus = 1/16;   % RasB
k5plus = 1/11;   % MHCKA

k1minus = 1/120;    % GBgamma  
k2minus = 1/130;    % GBPC
k3minus = 1/145;    % MCOR
k4minus = 1/160;    % RasB
k5minus = 1/115;    % MHCKA

savefigure = 0; % Set 1 if want to save figure and set 0 if do not want to save figure

% Put parameter values into vector

params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    k1minus, k2minus, k3minus, k4minus, k5minus];

% Set time span

tspan = [0 1200];

% Set initial conditions

init_cond = [0.5 0.5 0.5 0.5 0.5];

% Run Ode Solver

[T,y] = ode45(@(t,Y) blebSolverSteadyStates(t,Y,R,params) , tspan ,...
    init_cond);


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