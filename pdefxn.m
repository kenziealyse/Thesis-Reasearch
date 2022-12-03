function GBG_soln = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus)

% Function to solve PDE

% SET UP THE COMPUTATIONAL GRID
g = k1minus;
m = 9;           % number of interior nodes
numnodes = m+1;  % number of computational nodes
deltaX = 1/(m+1);     % uniform grid size in each cordinate
xnodes = linspace(0,lengthScale,numnodes+1);
alpha = d/deltaX^2;
time_step = final_time/deltaT;

% CREATE COEFFICIENT MATRIX 
e = ones(numnodes,1);
Amat = spdiags([alpha*e -(g+2*alpha)*e alpha*e],-1:1,numnodes,numnodes);
Amat(1,2) = 2*alpha; 
Bound = zeros(numnodes,1);
Bound(numnodes,1) = GBGSS*alpha;

k = lengthScale/2;
yint = GBGSS;
a = yint/(k^2);
l = 0;

init_solnmat = zeros(numnodes,1);
        
for i = 2:m+2
    init_solnmat(i-1,1) = a*(xnodes(i-1)-k)^2+l;    
end

% SOLVE THE ALGEBRAIC EQUATION
GBG_soln = zeros(1,time_step);

for i = 1:time_step
solution = (speye(numnodes) - deltaT.*Amat)\(deltaT.*Bound + init_solnmat);
init_solnmat = solution; 
GBG_soln(i) = solution(1);
end

% figure(3)
% plot(GBG_soln)

end

