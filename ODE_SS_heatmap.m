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
sigma_values = 6:.01:7.42;    %sigma


% Preallocate space
MHCKA_values = zeros(length(sigma_values), length(alpha_values));

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, ~, ~, ~, ~, ~,...
    k1minus, beta, ~, ~, ~, ~, k_0, ~]= setParameters();

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
        Y = SSsolnsSolver(R,params, k_0);

        GBG_0 = Y(1);
        GBPC_0 = Y(2);
        MCOR_0 = Y(3);
        RASB_0 = Y(4);
        MHCKA_0 = Y(5);

        MHCKA_values(p,i) = MHCKA_0;
    end 
end

min(MHCKA_values)

figure(1)
imagesc(alpha_values, sigma_values, MHCKA_values);
xlabel('\bf \alpha Value','FontSize',17);
ylabel('\bf \sigma Value','FontSize',17);
colorbar;
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/ODE_SS_smaller_range_','heatmap.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder