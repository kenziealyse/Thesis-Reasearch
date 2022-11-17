clear all; close all; clc

alpha1 = 50;
alpha2 = 15;

f = @(t,Y) [1 - Y(1)*(1 + Y(2)); alpha1 - alpha2*Y(2)];


y1 = linspace(0,100,20);
y2 = linspace(0,100,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 column

[x,y] = meshgrid(y1,y2);
size(x)
size(y)

u = zeros(size(x));
v = zeros(size(x));


t=0; % we want the derivatives at each point at t=0, i.e. the starting time


for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end


quiver(x,y,u,v,'r');
figure(gcf)
xlabel('X_{CR}')
ylabel('X_{MK}')
axis tight equal;
