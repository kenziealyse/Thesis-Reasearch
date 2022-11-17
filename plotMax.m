function plotMax(T, max_MhcA_A, max_MhcA_B, max_MhcA_C, str)
% Function to plot the max MhcA levels

plot(T, max_MhcA_A, 'linewidth', 3)
hold on
plot(T, max_MhcA_B, 'linewidth', 3)
hold on
plot(T, max_MhcA_C, 'linewidth', 3)

legend('Boundary A', 'Boundary B', 'Boundary C');
set(0,'defaultaxesfontsize',40);
xlabel('\bf Time (Secs)')
ylabel('\bf Relative Fluoresence')

end
