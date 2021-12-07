%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to Run Ode Solver using the blebSolver
%Declares variables, solves ODE
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