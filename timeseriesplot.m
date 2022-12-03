%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to Run Ode Solver and save a times series plot as a figure
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
close all
clear

savefigure = 1; % Set 1 if want to save figure and set 0 if do not want to save figure

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0] = setParameters();

% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    k1minus, k2minus, k3minus, k4minus, k5minus, k_0];

% Set time span
tspan = [0 200];

% Set initial conditions
init_cond = [.5 .01 .01 .01 .01];

% Run Ode Solver
[T,y] = ode45(@(t,Y) blebSolver(t,Y,R,params) , tspan ,...
    init_cond);

% Plot solutions
figure(1)
GBG = y(:,1);
GBPC = y(:,2);
MCOR = y(:,3);
RASB = y(:,4);
MHCKA = y(:,5); 
plot(T, GBG, 'LineWidth', 1, 'Color', 'r')
hold on
plot(T, GBPC, 'LineWidth', 2.5, 'Color', 'b', 'LineStyle','-.')
plot(T, MCOR, 'LineWidth', 2, 'Color', 'm', 'LineStyle',':')
plot(T, RASB, 'LineWidth', 1.5, 'Color', 'g', 'LineStyle','-')
plot(T, MHCKA, 'LineWidth', 1, 'Color', 'k', 'LineStyle','--')
set(0,'defaultaxesfontsize',16);
ylim([0 1])
xlabel("\bf Time (Seconds)")
ylabel("\bf Concentration")
legend('GBG', 'GBPC', 'MCOR', 'RasB', 'MHCKA','location','best')

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/SSSolnsPlot', '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder

for i = 1:5 
    SSsolns(i, 1) = y(end,i); 
end

VarNames = {'GBG Steady State', 'GBPC Steady State', 'MCOR Steady State',...
        'RASB Steady State', 'MHCKA Steady State'};

steady_state_table = table(SSsolns(1,1),SSsolns(2,1),SSsolns(3,1)...
    ,SSsolns(4,1),SSsolns(5,1), 'VariableNames',VarNames)


