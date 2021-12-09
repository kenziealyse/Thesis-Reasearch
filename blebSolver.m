function  dY_dt = blebSolver(~,Y,R,params)

%BLEBSOLVER is a function that holds all the equations for the differential
%equation and takes in the variables needed to solve the differential
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

dGBG_dt = alpha1*R*(1-GBG)- k1*GBG;          %eq (G beta gamma)
dGBPC_dt =  alpha2*GBG*(1-GBPC)-k2*GBPC;  %eq (GBPC)
dMCOR_dt =  alpha3*GBPC*(1-MCOR)-k3*MCOR*MHCKA;  %eq (XMCOR)
dRASB_dt =  alpha4*GBG*(1-RASB)-k4*RASB;  %eq (RASB)
dMHCKA_dt = alpha5*RASB*(1-MHCKA)-k5*MHCKA;  %eq (XMCHKA)

dY_dt = [dGBG_dt; dGBPC_dt; dMCOR_dt; dRASB_dt; dMHCKA_dt];

end