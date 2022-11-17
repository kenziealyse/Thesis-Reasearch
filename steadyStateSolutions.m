function [T,SSsolns] = steadyStateSolutions(R, params)
% Function to Run Ode Solver and look at the steady states
% when changing variables
% Varaible changed: R0
% Returned: SS Solutions, Time

% Set time span
tspan = [0 200];

% Set initial conditions
init_cond = [.5 0 0 0 0];

% Pre-allocate Space in SSsolns vector
SSsolns = zeros(length(R), length(init_cond));

% Run Ode Solver in loop for various R0 values
for i = 1:length(R)
    
    [T,y] = ode45(@(t,Y) blebSolver(t,Y,R(i),params) , tspan ,...
    init_cond);
    
    SSsolns(i, :) = y(end, :); 
end

