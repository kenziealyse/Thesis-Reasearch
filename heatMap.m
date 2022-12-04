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

% Set alpha and beta values for for loop
alpha_values = 5:.5:60;    %alpha
sigma_values = .5:.01:7;    %sigma


% Preallocate space
percentChange = zeros(length(sigma_values), length(alpha_values));

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, ~, ~, ~, ~, ~,...
    k1minus, beta, ~, ~, ~, ~, k_0] = setParameters();

for i = 1:length(alpha_values)
    alpha = alpha_values(i);         % alpha
    k2plus = alpha*k1plus;           % GBPC
    k3plus = alpha*k1plus;           % MCOR
    k4plus = alpha*k1plus;           % RasB
    k5plus = alpha*k1plus;           % MHCKA
    k3minus = alpha*k1plus;          % MCOR
    
    for p = 1:length(sigma_values)
        beta = sigma_values(p)*alpha_values(i);                  % beta
        k2minus = beta*k1minus;          % GBPC
        k4minus = beta*k1minus;          % RasB
        k5minus = beta*k1minus;          % MHCKA
        k6minus = beta*k1minus;          % MCOR degradation
        
        % Param Names
        params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
        k1minus, k2minus, k3minus, k4minus, k5minus];
        % Set time span
        tspan = 0:deltaT:final_time;
        Y = blebSolverforPDE(R,params, final_time, d, lengthScale, deltaT, ...
        tspan, k_0);

        GBPCprime = Y(:,1);
        MCORprime = Y(:,2);
        RasBprime = Y(:,3);
        MHCKAprime = Y(:,4);
        % Calculate the percent change of myosin
        percentChange(p,i) = 100*((max(MCORprime) - MCORprime(1))/MCORprime(1));

    end 
end

figure(1)
imagesc(alpha_values, sigma_values, percentChange);
xlabel('\bf \alpha Value','FontSize',17);
ylabel('\bf \sigma Value','FontSize',17);
colorbar;
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
hold on
plot(30,5.5,'.', 'Color', 'r', 'MarkerSize', 20)

rectangle('Position', [15, 4.5, 20, 2], 'EdgeColor','r', 'LineWidth', 1.5) 

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/parameter','heatmap_with_box_and_point.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder