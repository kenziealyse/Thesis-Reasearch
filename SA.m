%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file plots a histogram when
% varying alpha or beta (specified by
% a string) for all activation and
% deactivation rates, respecively,
% within the specified range.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
close all
clear
clc

% Which value do you want to randomize? (alpha or beta)
str = 'alpha';
num_of_rand_sample = 10000;

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
    index = index_values(i); 
    var1 = var1_values(i); 
    % Figure Names 
    figureName = ["k2plus", "k3plus", "k4plus", "k5plus", "k2minus", ...
        "k3minus", "k4minus", "k5minus"];
        
    % Set paramets values
    [R, lengthScale, d, deltaT, final_time,...
        k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
        k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
        k_0, ~] = setParameters();
    
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
    end
    
    % Specifiy histogram limits and edges
    myxlim = [-30,30];
    edges = -30:1:30;

    % Plot the histogram
    figure(i)
    histogram(percentChange, edges);
    xlabel('\bf Percent Change','FontSize',17);
    ylabel('\bf Frequency','FontSize',17);
    xlim(myxlim)
    xticks(-30:5:30)
    
    % Save the figure
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
    figure_name = ['/SA', figureName{index}, '_' str, '_randomized', 'histo.pdf'];
    dirPath = strcat('/','figures', figure_name); % Directory Path
    saveas(gcf,[pwd dirPath]); % Save Figure in Folder
end