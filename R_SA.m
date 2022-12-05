%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mackenzie Dalton
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
close all
clear
clc

num_of_rand_sample = 10000;

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0, ~] = setParameters();

% Set paramets values
init_R_testing = R;
lower_bound = 0.1;
upper_bound = 1;
R_values = lower_bound:.0000001:upper_bound;  
R_rand_sample = randsample(R_values, num_of_rand_sample);

% Set time span
tspan = 0:deltaT:final_time;

% Param Names
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    k1minus, k2minus, k3minus, k4minus, k5minus];

% Find initial percent change
Y = blebSolverforPDE(R,params, final_time, d, lengthScale, deltaT, ...
    tspan, k_0);
GBPCprime = Y(:,1);
MCORprime = Y(:,2);
RasBprime = Y(:,3);
MHCKAprime = Y(:,4);
initial_percent = 100*((max(MCORprime) - MCORprime(1))/MCORprime(1));

% Pre Allocate Space 
percentChange = zeros(length(R_rand_sample),1);

for jj = 1:length(R_rand_sample)
    R = R_rand_sample(jj);           
    Y = blebSolverforPDE(R, params, final_time, d, lengthScale, deltaT, ...
    tspan, k_0);
    GBPCprime = Y(:,1);
    MCORprime = Y(:,2);
    RasBprime = Y(:,3);
    MHCKAprime = Y(:,4);
    % Calculate the percent change of myosin
    percentChange(jj) = initial_percent - 100*((max(MCORprime) - MCORprime(1))/MCORprime(1));
end


myxlim = [-30,30];
edges = -30:1:30;

figure(1)
histogram(percentChange, edges);
xlabel('\bf Percent Change','FontSize',17);
ylabel('\bf Frequency','FontSize',17);
xlim(myxlim)
xticks(-30:5:30);

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = '/SA_R_randomized_histo.pdf';
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder