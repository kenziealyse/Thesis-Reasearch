function  dY_dt = blebSolverSteadyStates(t,Y,k5plusval,params, GBG, GBPC, RASB)

%BLEBSOLVERSTEADYSTATES is a function that holds all the equations for the differential
%equation at the steady state solutions and takes in the variables needed to solve the differential
%equations


% Set k plus paramets

k3plus = params(1);


% Set K minus parameters

k3minus = params(2);
k5minus = params(3);


% Relabel to easily keep track of compartments

MCOR = Y(1);
MHCKA = Y(2);

k5plus = k5plusval(t);       

dMCOR_dt =  (k3plus*GBG)/(k3plus*GBPC+k3minus*MHCKA);  % Steady State eq (XMCOR)
dMHCKA_dt = (k5plus*RASB)/(k5plus*RASB-k5minus);  % Steady State eq (XMCHKA)

dY_dt = [dMCOR_dt; dMHCKA_dt];

end

