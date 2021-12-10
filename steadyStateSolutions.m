%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to Run Ode Solver and look at the steady states
% when changing variables
% Varaible changed: R0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE

close all
clear all

% Set paramets values

R = (0:.1:1);          % Initial active receptors 
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

savefigure = 1; % Set 1 if want to save figure and set 0 if do not want to save figure

% Put parameter values into vector

params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    k1minus, k2minus, k3minus, k4minus, k5minus];

% Set time span

tspan = [0 120];

% Set initial conditions

init_cond = [.5 0 0 0 0];

% Pre-allocate Space in SSsolns vector

SSsolns = zeros(length(R), length(init_cond));

% Run Ode Solver in loop for various R0 values

for i = 1:length(R)
    
    [T,y] = ode45(@(t,Y) blebSolver(t,Y,R(i),params) , tspan ,...
    init_cond);

    % Save Steady State Solutions 

    
    SSsolns(i, :) = y(end); 


end

T(end)

VarNames = {'R_0', 'GBG Steady State', 'GBPC Steady State', 'MCOR Steady State',...
        'RASB Steady State', 'MHCKA Steady State'};

steady_state_table = table(R', SSsolns(:,1),SSsolns(:,2),SSsolns(:,3)...
    ,SSsolns(:,4),SSsolns(:,5), 'VariableNames',VarNames)


