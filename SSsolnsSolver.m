function  SSsolns = SSsolnsSolver(R,params, k_0)

%SSsolnsSolver is a function that holds all the S.S. equations for the 
%differential equations and takes in the parameters needed to solve 
%the differential equations and returns the S.S. solutions

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

% System of equations
GBG_0 =   (k1plus*R)/(k1plus*R + k1minus);                    %eq (G beta gamma)
GBPC_0 =  (k2plus*GBG_0)/(k2plus*GBG_0 + k2minus);            %eq (GBPC)
RASB_0 =  (k4plus*GBG_0)/(k4plus*GBG_0 + k4minus);            %eq (RASB)
MHCKA_0 = 1 - (k5minus/(k5plus*(RASB_0 + k_0)));              %eq (XMCHKA)
MCOR_0 =  (k3plus*GBPC_0)/(k3plus*GBPC_0 + k3minus*MHCKA_0);  %eq (XMCOR)

% Solutions vector to return
SSsolns = [GBG_0; GBPC_0; MCOR_0; RASB_0; MHCKA_0];

end