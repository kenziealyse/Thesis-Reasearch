function max_sigma = sigma(plotn)

%sigma is a function that calculates the max sigma value based on the given
%parameter values. It also has a plotn input that is 1 to plot he desired 
%region for beta\alpha and 0 otherwise.

% Set paramets values
[R, ~, ~, ~, ~, k1plus, ~, ~, ~, ~, ~, ~,k1minus, ~, ~, ~, ~, ~,~]...
    = setParameters();

% Set gamma and phi
gamma = (k1minus/k1plus)*(1 + (k1minus/k1plus)*(1/R));
phi = k1plus/k1minus;

% Calculate the max sigma value
max_sigma = ((-1 + sqrt(1 + 4*gamma*phi)))/(2*gamma);

% Plot the desired region for beta\alpha
if plotn == 1
x = -20:.001:10;
x_vals = 0:.001:max_sigma;
figure(1)
plot(x, gamma.*x.^2 + x - phi, 'LineWidth', 1.5)
hold on
plot(x_vals, gamma.*x_vals.^2 + x_vals - phi, 'LineWidth', 1.5, 'Color', 'r')
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
end

end

