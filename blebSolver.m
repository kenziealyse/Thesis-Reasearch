function  dY_dt = blebSolver(~,Y,R,params)

%BLEBSOLVER is a function that holds all the equations for the differential
%equation and takes in the variables needed to solve the differential
%equations
%testing

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
k_0 = params(11);

% Relabel to easily keep track of compartments
GBG = Y(1);
GBPC = Y(2);
MCOR = Y(3);
RASB = Y(4);
MHCKA = Y(5); 

% System of equations
dGBG_dt =   k1plus*R*(1-GBG) - k1minus*GBG;                            %eq (G beta gamma)
dGBPC_dt =  k2plus*GBG*(1-GBPC)- k2minus*GBPC;  %eq (GBPC)
dMCOR_dt =  k3plus*GBPC*(1-MCOR)- k3minus*MCOR*MHCKA;  %eq (XMCOR)
dRASB_dt =  k4plus*GBG*(1-RASB) - k4minus*RASB;  %eq (RASB)
dMHCKA_dt = k5plus*(RASB+k_0)*MHCKA*(1-MHCKA) - k5minus*MHCKA;  %eq (XMCHKA)

% Solutions vector to return
dY_dt = [dGBG_dt; dGBPC_dt; dMCOR_dt; dRASB_dt; dMHCKA_dt];

end