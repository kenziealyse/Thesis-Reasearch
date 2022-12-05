%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0, myxlim] = setParameters();

% Set time span
tspan = 0:deltaT:final_time;

% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    k1minus, k2minus, k3minus, k4minus, k5minus];

for i = 1:length(R)  
Y = blebSolverforPDE(R(i),params, final_time, d, lengthScale, deltaT, ...
    tspan, k_0);
end

GBPCprime = Y(:,1);
MCORprime = Y(:,2);
RasBprime = Y(:,3);
MHCKAprime = Y(:,4);

% Calculate the percent change of myosin
percentChange = 100*((max(MCORprime) - MCORprime(1))/MCORprime(1))

% Save the figures
myosin_fig = figure(2);
plot(tspan, MCORprime, 'LineWidth', 2, 'Color', 'm', 'LineStyle',':')
xlim(myxlim)
xlabel('\bf Time (Seconds)', 'Fontsize', 17)
ylabel('\bf Concentration', 'Fontsize', 17)
set(myosin_fig, 'Units', 'Inches');
pos = get(myosin_fig, 'Position');
set(myosin_fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
figure_name = ['/MyosinTimeSeriesPlot_redpoint', '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder

all_fig = figure(1);
plot(tspan, GBPCprime, 'LineWidth', 2.5, 'Color', 'b', 'LineStyle','-.')
hold on
plot(tspan, MCORprime, 'LineWidth', 2, 'Color', 'm', 'LineStyle',':')
plot(tspan, RasBprime, 'LineWidth', 1.5, 'Color', 'g', 'LineStyle','-')
plot(tspan, MHCKAprime, 'LineWidth', 1, 'Color', 'k', 'LineStyle','--')
legend('GBPC','Myosin', 'RasB', 'MHCKA', 'Location', 'Best', 'FontSize', 17)
ylabel('\bf Concentration', 'FontSize', 17)
xlabel('\bf Time (Seconds)', 'FontSize', 17)
xlim(myxlim)
set(all_fig, 'Units', 'Inches');
pos = get(all_fig, 'Position');
set(all_fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
figure_name = ['/ALLTimeSeriesPlot_redpoint', '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder


