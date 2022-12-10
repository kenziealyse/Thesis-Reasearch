function  Y = blebSolverforPDE(R,params, final_time, d, lengthScale, deltaT, ...
    tspan, k_0)

%BLEBSOLVERforPDE is a function that solves the differential equations with
%the pde solution. It first finds the steadt state soltion of each species,
%then it solves the pde, and lastly it uses the output of the pde to solve
%the ODE system. To solve the ODE system, a forward difference scheme is 
%used. It returns the solutions to the ODE system in vector
%form.

    % Preallocate Space
    RasBprime = zeros(length(tspan), 1);
    GBPCprime = zeros(length(tspan), 1);
    MCORprime = zeros(length(tspan), 1);
    MHCKAprime = zeros(length(tspan), 1);

    % Label parameter values
    k2plus = params(2);
    k3plus = params(3);
    k4plus = params(4);
    k5plus = params(5);
    k1minus = params(6);
    k2minus = params(7);
    k3minus = params(8);
    k4minus = params(9);
    k5minus = params(10);

    % Calculate Steady State Values
    SSsolns = SSsolnsSolver(R, params, k_0);
    GBGSS =  SSsolns(1);  % Steady State Value of GBG

    % Use PDE for GBG
    GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);

    % Set initial conditions
    init_cond = [SSsolns(2) SSsolns(3) SSsolns(4)...
       SSsolns(5)];
    GBPCprime(1) = init_cond(1);
    MCORprime(1) = init_cond(2);
    RasBprime(1) = init_cond(3);
    MHCKAprime(1) = init_cond(4);

    for jj = 1:length(tspan) - 1   
         %eq (GBPC)
         GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) ...
             - k2minus*GBPCprime(jj))*deltaT;
         %eq (XMCOR)
         MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
            k3minus*MCORprime(jj)*MHCKAprime(jj))*deltaT;
         % eq (RASB)
         RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
             k4minus*RasBprime(jj))*deltaT; 
         %eq (XMCHKA)
         MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*MHCKAprime(jj)*...
            (RasBprime(jj)+ k_0)*(1 - MHCKAprime(jj))...
              - k5minus*MHCKAprime(jj))*deltaT;  
    end  
    
    % Solutions Vector
    Y = [GBPCprime MCORprime RasBprime MHCKAprime];
end

