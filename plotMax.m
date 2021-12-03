function plotMax(T, max_MhcA_A, max_MhcA_B, max_MhcA_C, str)
 

plot(T, max_MhcA_A, 'linewidth', 3)
hold on
plot(T, max_MhcA_B, 'linewidth', 3)
hold on
plot(T, max_MhcA_C, 'linewidth', 3)

legend('Boundary A', 'Boundary B', 'Boundary C');
ylabel('Fluoresence (AU)')
xlabel('Time')
title('Max Myosin Levels Over Time for', str)


hello

chqnge

end
