function  dY_dt = blebSolverSteadyStates(~,Y,R,params)

%BLEBSOLVERSTEADYSTATES is a function that holds all the equations for the differential
%equation at the steady state solutions and takes in the variables needed to solve the differential
%equations


% Set k plus paramets

k1plus = params(1);
k2plus = params(2);
k3plus = params(3);
k4plus = params(4);
k5plus = params(5);


% Set K minus parameters

k1minus = params(6);
k2minus = params(7);
k3minus = params(8);
k4minus = params(9);
k5minus = params(10);


% Relabel to easily keep track of compartments
GBG = Y(1);
GBPC = Y(2);
MCOR = Y(3);
RASB = Y(4);
MHCKA = Y(5);

dGBG_dt = (k1plus*R)/(k1plus*R+k1minus);          %eq (G beta gamma)
dGBPC_dt =  (k2plus*GBG)/(k2plus*GBG+k2minus);  %eq (GBPC)
dMCOR_dt =  (k3plus*GBG)/(k3plus*GBPC+k3minus*MHCKA);  %eq (XMCOR)
dRASB_dt =  (k4plus*GBG)/(k4plus*GBG-k4minus);  %eq (RASB)
dMHCKA_dt = (k5plus*RASB)/(k5plus*RASB-k5minus);  %eq (XMCHKA)

dY_dt = [dGBG_dt; dGBPC_dt; dMCOR_dt; dRASB_dt; dMHCKA_dt];

end

