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

% Which value do you want to randomize? (alpha or beta)
str = 'alpha';
num_of_rand_sample = 5001;

if strcmp(str,'beta')
var1_values = [7, 8, 9, 10];
index_values = [5, 6, 7, 8];
random_index = 1;
elseif strcmp(str,'alpha')
var1_values = [2, 3, 4, 5];
index_values = [1, 2, 3, 4];
random_index = 2;
end

for i = 1:length(var1_values)
% SET THESE

index = index_values(i); % which string do you want?
var1 = var1_values(i); % which variable value in the legend? [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    %k1minus, k2minus, k3minus, k4minus, k5minus];

% Figure Names 
figureName = ["k2plus", "k3plus", "k4plus", "k5plus", "k2minus", ...
    "k3minus", "k4minus", "k5minus"];

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0, myxlim] = setParameters();

randomized_params = [beta, alpha];
randomized_params2 = [k1minus, k1plus];

% Set paramets values
init_var_testing = randomized_params(random_index);
lower_bound = init_var_testing - .2*init_var_testing;
upper_bound = init_var_testing + .2*init_var_testing;
var_values = lower_bound:.00001:upper_bound;  % alpha
var_rand_sample = randsample(var_values, num_of_rand_sample);

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
percentChange = zeros(length(var_rand_sample),1);
MCORprime2 = zeros(length(tspan),length(var_rand_sample));

for jj = 1:length(var_rand_sample)
    params(var1) = var_rand_sample(jj)*randomized_params2(random_index);           
    Y = blebSolverforPDE(R, params, final_time, d, lengthScale, deltaT, ...
    tspan, k_0);
    GBPCprime = Y(:,1);
    MCORprime = Y(:,2);
    RasBprime = Y(:,3);
    MHCKAprime = Y(:,4);
    % Calculate the percent change of myosin
    percentChange(jj) = initial_percent - 100*((max(MCORprime) - MCORprime(1))/MCORprime(1));
    MCORprime2(:,jj) = MCORprime;
end

figure(i)
[lineh, bandsh] = fanChart(tspan, MCORprime2, 'median');
txt = strcat({'Pct'}, cellstr(int2str((10:10:90)')));
legend([lineh;bandsh], [{'Median'};txt], 'location', 'best')
ylabel('\bf Concentraion','FontSize',17)
xlabel('\bf Time (seconds)','FontSize',17)
xlim(myxlim)
ylim([.5 1])
yticks(.5:.1:1)

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/PrctilePlot', figureName{index}, '_' str, '.pdf'];
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder
end