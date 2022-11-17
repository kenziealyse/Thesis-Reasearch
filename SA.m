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

% SET THESE

index = 1; % which string do you want?
var1 = 8; % which variable value in the legend? params = [k1plus, k2plus,
% k3plus, k4plus, k5plus, k1minus, k2minus, k3minus, k4minus, k5minus];

% Figure Names 
figureName = ['k2plus', 'k3plus', 'k4plus', 'k5plus', 'k3minus'];

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0] = setParameters();

% Set paramets values
init_alpha = 2;
lowerBound = init_alpha - .3*init_alpha;
upperBound = init_alpha + .3*init_alpha;
alpha_values = lowerBound:.001:upperBound;  % alpha
alpha_rand_sample = randsample(alpha_values, 1000);
% myxlim = [-12, 12];
% edges = -12:1:12;

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


for i = 1:length(alpha_rand_sample)
    alpha = alpha_rand_sample(i);
    params(var1) = alpha*params(var1);           % GBPC
    Y = blebSolverforPDE(R, params, final_time, d, lengthScale, deltaT, ...
    tspan, k_0);
    GBPCprime = Y(:,1);
    MCORprime = Y(:,2);
    RasBprime = Y(:,3);
    MHCKAprime = Y(:,4);

    % Calculate the percent change of myosin
    percentChange(i) = initial_percent - 100*((max(MCORprime) - MCORprime(1))/MCORprime(1));
end

figure(1)
histogram(percentChange);
xlabel('\bf Percent Change','FontSize',17);
ylabel('\bf Frequency','FontSize',17);
% xlim(myxlim)

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/SA', figureName(index), 'histo.pdf']  
dirPath = strcat('/','figures', figure_name); % Directory Path
%saveas(gcf,[pwd dirPath]); % Save Figure in Folder