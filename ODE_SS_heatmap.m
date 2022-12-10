%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file plots a heat map of the 
% numerical and calculated values of
% the steady state solution for MHCKA.
% The heat map is broken up into 
% positive, 0, and negative.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
close all
clear
clc

% Set alpha and beta values for for loop
alpha_values = 5:1:15;         % alpha
sigma_values = .5:.0001:20;    % sigma
m = 2;                         % 10^(m) tolerance
% Set initial conditions
init_cond = [.5 .01 .01 .01 0.1];
% Set time span
tspan = [0 200];

% Preallocate space
MHCKA_values_numerical = zeros(length(sigma_values), length(alpha_values));
MHCKA_values_calculated = zeros(length(sigma_values), length(alpha_values));

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
        k1minus, k2minus, k3minus, k4minus, k5minus, k_0];
        % Run Ode Solver
        [T,y] = ode45(@(t,Y) blebSolver(t,Y,R,params) , tspan ,...
        init_cond);
        SSsolns = SSsolnsSolver(R,params, k_0);


        MHCKA_numerical_0 = y(end,5);
        MHCKA_calculated_0 = SSsolns(5);
        
       % Check how close steady state is to 0 
       if MHCKA_calculated_0 > 10^(-m) && MHCKA_calculated_0 < 0.1
                MHCKA_calculated_0 = 0.1;
       end
       if MHCKA_calculated_0 < -10^(-m) && MHCKA_calculated_0 > -0.1
                MHCKA_calculated_0 = -0.1;
       end

       if MHCKA_numerical_0 > 10^(-m) && MHCKA_numerical_0 < 0.1
                MHCKA_numerical_0 = 0.1;
       end
       if MHCKA_numerical_0 < -10^(-m) && MHCKA_numerical_0 > -0.1
                MHCKA_numerical_0 = -0.1;
       end

       % Save steady state values
       MHCKA_values_numerical(p,i) = MHCKA_numerical_0;
       MHCKA_values_calculated(p,i) = MHCKA_calculated_0;
    end 
end

% Set color map, split colors by negative, 0 and positive
cmap = [247 239 210
        247 239 210
        247 239 210
        247 239 210
        71 92 108
        205 139 98
        205 139 98
        205 139 98
        205 139 98];

% switch colors to RGB numbers
cmap = cmap./255;

% Plot the heat map for numerical values
figure(1)
imagesc(alpha_values, sigma_values, MHCKA_values_numerical);
hold on
yline(7.2129, 'linewidth', 2, 'Color', 'm')
hold off
xlabel('\bf \alpha Value','FontSize',17);
ylabel('\bf \sigma Value','FontSize',17);
colormap(cmap);
colorbar('XTickLabel',{'-0.4','-0.3','-0.2','-0.1',...
               '0','0.1','0.2','0.3','0.4'}, ...
               'XTick', -.4:.1:.4)
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
set(gcf, 'Units', 'Inches');
set(gca, 'CLim', [-.4, .4]);

% Save the figure
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
figure_name = ['/ODE_SS_numerical','heatmap.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder

% Plot the heat map for calculated values values
figure(2)
imagesc(alpha_values, sigma_values, MHCKA_values_calculated);
colormap(parula(5))
hold on
yline(7.2129, 'linewidth', 2.5, 'Color', 'm')
xlabel('\bf \alpha Value','FontSize',17);
ylabel('\bf \sigma Value','FontSize',17);
colormap(cmap);
set(gca, 'CLim', [-.4, .4]);
colorbar('XTickLabel',{'-0.4','-0.3','-0.2','-0.1',...
               '0','0.1','0.2','0.3','0.4'}, ...
               'XTick', -.4:.1:.4)
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented

% Save the figure
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
figure_name = ['/ODE_SS_calculated','heatmap.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder