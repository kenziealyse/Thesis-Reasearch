clc
clear 
close all

[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0, myxlim] = setParameters();
% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
    k1minus, k2minus, k3minus, k4minus, k5minus, k_0];


SSsolns = SSsolnsSolver(R,params, k_0);
RasSS = SSsolns(4);

bifpoint = 7.48;

x1 = 0:.01:bifpoint;
y1 = 1 - x1.*k1minus/(k1plus*(RasSS +k_0));
z1 = 0.*x1;

x2 = bifpoint:.01:30;
y2 = 1 - x2.*k1minus/(k1plus*(RasSS +k_0));
z2 = 0.*x2;

figure(1)
plot(x1,y1, 'linewidth', 2.5, 'Color', [205 139 98]./255)
hold on
plot(x2,y2, 'LineWidth', 2.5, 'Color', [238 215 161]./255, 'LineStyle','-.')
plot(x1,z1, 'LineWidth', 2.5, 'Color', [71 92 108]./255, 'LineStyle','-.')
plot(x2,z2, 'linewidth', 2.5, 'Color', [71 92 108]./255)
grid on
xline(7.2129,'linewidth', 1.5, 'color', 'm')
xlim([0 15])
ylim([-1 1])
xlabel('\bf \sigma Value','FontSize',17);
ylabel('\bf MHCKA S.S. Value','FontSize',17);
legend('Stable non-trivial SS', 'Unstable non-trivial SS', ...
    'Unstable trivial SS', 'Stable trivial SS', 'Sufficient Positivity Criteria')


set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/BifurcationPlot', '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder