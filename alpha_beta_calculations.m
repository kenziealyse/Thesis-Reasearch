% Clear the workspace
clc
clear
close all

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0] = setParameters();

gamma = (k1minus/k1plus)*(1 + (k1minus/k1plus)*(1/R));
sigma = k1plus/k1minus;

max_beta = ((-1 + sqrt(1 + 4*gamma*sigma)))/(2*gamma)

x = -20:.001:10;
x_vals = 0:.001:max_beta;

figure(1)
plot(x, gamma.*x.^2 + x - sigma, 'LineWidth', 1.5)
hold on
plot(x_vals, gamma.*x_vals.^2 + x_vals - sigma, 'LineWidth', 1.5, 'Color', 'r')
yline(0, 'LineWidth', 1.5)
xline(0, 'LineWidth', 1.5)
xlabel('\bf $\frac{\beta}{\alpha}$','FontSize',25, 'Interpreter','latex');
ylabel('y','FontSize',20);
grid on


set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
figure_name = '/positivity_plot.pdf';   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder

