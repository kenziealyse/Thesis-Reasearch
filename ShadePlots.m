%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mackenzie Dalton
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
close all
clear
clc

%% Set Global Params

% Set alpha and beta values for for loop
R = 0.8;                         % Initial active receptors 
lengthScale = 2.725;             % Diffusion length scale
d = .011;                        % Diffusion Rate
deltaT = .01;                    % Time Step
final_time = 13;                 % Final Time
init_alpha = 15;
beta = 5;
lowerBound = init_alpha - .3*init_alpha;
upperBound = init_alpha + .3*init_alpha;
dlowerBound = d - .3*d;
dupperBound = d + .3*d;
d_values = dlowerBound:.00001:dupperBound;
alpha_values = lowerBound:.001:upperBound;  % alpha
alpha_rand_sample = randsample(alpha_values, 500);
R_rand_sample = randsample(0:.0001:1, 500);
d_rand_sample = randsample(d_values, 50);
k1plus = 1/10;                           % GBgamma
k1minus = 1/120;                         % GBgamma  
k2minus = beta*k1minus;          % GBPC
k4minus = beta*k1minus;          % RasB
k5minus = beta*k1minus;          % MHCKA
k6minus = beta*k1minus;          % MCOR degradation



%% k2plus

figureName = 'k2plus';

% Fixed Parameters
k3plus = init_alpha*k1plus;           % MCOR
k4plus = init_alpha*k1plus;           % RasB
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus;          % MCOR



for i = 1:length(alpha_rand_sample)
    alpha = alpha_rand_sample(i);
    k2plus = alpha*k1plus;           % GBPC
      
        % Set time span
        tspan = 0:deltaT:final_time;
        
        % Preallocate Space
        RasBprime = zeros(length(tspan), 1);
        GBPCprime = zeros(length(tspan), 1);
        MCORprime = zeros(length(tspan), 1);
        MHCKAprime = zeros(length(tspan), 1);
        
        % Put parameter values into vector
        params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
        k1minus, k2minus, k3minus, k4minus, k5minus];
        
        % Calculate Steady State Values
        SSsolns = steadyStateSolutions(R, params);
        GBGSS =  SSsolns(:,1);  % Steady State Value of GBG
        
        % Use PDE for GBG
        GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);
        
        % Set initial conditions
        init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
           SSsolns(:,5)];
        GBPCprime(1) = init_cond(1);
        MCORprime(1) = init_cond(2);
        RasBprime(1) = init_cond(3);
        MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
        figure(1)
        plot(tspan, MCORprime, 'Color', [.7 .7 .7], 'linewidth', 2)
        hold on
end
      
% Set time span
tspan = 0:deltaT:final_time;

% Fixed Parameters
k2plus = init_alpha*k1plus;           % MCOR
k3plus = init_alpha*k1plus;           % MCOR
k4plus = init_alpha*k1plus;           % RasB
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus;          % MCOR

% Preallocate Space
RasBprime = zeros(length(tspan), 1);
GBPCprime = zeros(length(tspan), 1);
MCORprime = zeros(length(tspan), 1);
MHCKAprime = zeros(length(tspan), 1);

% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
k1minus, k2minus, k3minus, k4minus, k5minus];

% Calculate Steady State Values
SSsolns = steadyStateSolutions(R, params);
GBGSS =  SSsolns(:,1);  % Steady State Value of GBG

% Use PDE for GBG
GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);
        
% Set initial conditions
init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
   SSsolns(:,5)];
GBPCprime(1) = init_cond(1);
MCORprime(1) = init_cond(2);
RasBprime(1) = init_cond(3);
MHCKAprime(1) = init_cond(4);

for jj = 1:length(tspan)-1    
    GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
        k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
    MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
        k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
    RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
        k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
    MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
         - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
end

plot(tspan, MCORprime, 'Color', 'm', 'linewidth', 2)
xlabel('\bf Time (Seconds)', 'FontSize', 17)
ylabel('\bf Concentration', 'FontSize', 17)

% Save Plot

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/SAshadedPlots', figureName, '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder


%% k3minus SA

% Set figure Name
figureName = 'k3minus';

% Preallocate space
percentChange = zeros(1,length(alpha_rand_sample));

% Fixed Parameters
k2plus = init_alpha*k1plus;           % GBPC
k3plus = init_alpha*k1plus;           % RasB
k4plus = init_alpha*k1plus;           % MHCKA
k5plus = init_alpha*k1plus;           % MHCKA


for i = 1:length(alpha_rand_sample)
    alpha = alpha_rand_sample(i);
    k3minus = alpha*k1plus;           % GBPC
      
        % Set time span
        tspan = 0:deltaT:final_time;
        
        % Preallocate Space
        RasBprime = zeros(length(tspan), 1);
        GBPCprime = zeros(length(tspan), 1);
        MCORprime = zeros(length(tspan), 1);
        MHCKAprime = zeros(length(tspan), 1);
        
        % Put parameter values into vector
        params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
        k1minus, k2minus, k3minus, k4minus, k5minus];
        
        % Calculate Steady State Values
        SSsolns = steadyStateSolutions(R, params);
        GBGSS =  SSsolns(:,1);  % Steady State Value of GBG
        
        % Use PDE for GBG
        GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);
        
        % Set initial conditions
        init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
           SSsolns(:,5)];
        GBPCprime(1) = init_cond(1);
        MCORprime(1) = init_cond(2);
        RasBprime(1) = init_cond(3);
        MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
        % Calculate the percent change of myosin
        figure(2)
        plot(tspan, MCORprime, 'color', [.7 .7 .7], 'linewidth', 2)
        hold on
end


% Fixed Parameters
k2plus = init_alpha*k1plus;           % GBPC
k3plus = init_alpha*k1plus;           % RasB
k4plus = init_alpha*k1plus;           % MHCKA
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus; 
     
% Set time span
tspan = 0:deltaT:final_time;

% Preallocate Space
RasBprime = zeros(length(tspan), 1);
GBPCprime = zeros(length(tspan), 1);
MCORprime = zeros(length(tspan), 1);
MHCKAprime = zeros(length(tspan), 1);

% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
k1minus, k2minus, k3minus, k4minus, k5minus];

% Calculate Steady State Values
SSsolns = steadyStateSolutions(R, params);
GBGSS =  SSsolns(:,1);  % Steady State Value of GBG

% Use PDE for GBG
GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);

% Set initial conditions
init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
   SSsolns(:,5)];
GBPCprime(1) = init_cond(1);
MCORprime(1) = init_cond(2);
RasBprime(1) = init_cond(3);
MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
% Calculate the percent change of myosin
plot(tspan, MCORprime, 'color', 'm', 'linewidth', 2);
xlabel('\bf Time (Seconds)', 'FontSize', 17)
ylabel('\bf Concentration', 'FontSize', 17)

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/SAshadedPlots', figureName, '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder


%% R SA

% Set figure Name
figureName = 'R';

% Fixed Parameters
k2plus = init_alpha*k1plus;           % GBPC
k3plus = init_alpha*k1plus;           % RasB
k4plus = init_alpha*k1plus;           % MHCKA
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus;          % GBPC

for i = 1:length(R_rand_sample)
      
        % Set time span
        tspan = 0:deltaT:final_time;
        
        % Preallocate Space
        RasBprime = zeros(length(tspan), 1);
        GBPCprime = zeros(length(tspan), 1);
        MCORprime = zeros(length(tspan), 1);
        MHCKAprime = zeros(length(tspan), 1);
        
        % Put parameter values into vector
        params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
        k1minus, k2minus, k3minus, k4minus, k5minus];
        
        % Calculate Steady State Values
        SSsolns = steadyStateSolutions(R_rand_sample(i), params);
        GBGSS =  SSsolns(:,1);  % Steady State Value of GBG
        
        % Use PDE for GBG
        GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);
        
        % Set initial conditions
        init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
           SSsolns(:,5)];
        GBPCprime(1) = init_cond(1);
        MCORprime(1) = init_cond(2);
        RasBprime(1) = init_cond(3);
        MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
        % Calculate the percent change of myosin
        figure(3)
        plot(tspan, MCORprime, 'color', [.7 .7 .7], 'linewidth', 2)
        hold on
end


% Fixed Parameters
k2plus = init_alpha*k1plus;           % GBPC
k3plus = init_alpha*k1plus;           % RasB
k4plus = init_alpha*k1plus;           % MHCKA
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus; 
     
% Set time span
tspan = 0:deltaT:final_time;

% Preallocate Space
RasBprime = zeros(length(tspan), 1);
GBPCprime = zeros(length(tspan), 1);
MCORprime = zeros(length(tspan), 1);
MHCKAprime = zeros(length(tspan), 1);

% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
k1minus, k2minus, k3minus, k4minus, k5minus];

% Calculate Steady State Values
SSsolns = steadyStateSolutions(R, params);
GBGSS =  SSsolns(:,1);  % Steady State Value of GBG

% Use PDE for GBG
GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);

% Set initial conditions
init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
   SSsolns(:,5)];
GBPCprime(1) = init_cond(1);
MCORprime(1) = init_cond(2);
RasBprime(1) = init_cond(3);
MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
% Calculate the percent change of myosin
plot(tspan, MCORprime, 'color', 'm', 'linewidth', 2);
xlabel('\bf Time (Seconds)', 'FontSize', 17)
ylabel('\bf Concentration', 'FontSize', 17)

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/SAshadedPlots', figureName, '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder



%% d SA

% Set figure Name
figureName = 'd';

% Fixed Parameters
k2plus = init_alpha*k1plus;           % GBPC
k3plus = init_alpha*k1plus;           % RasB
k4plus = init_alpha*k1plus;           % MHCKA
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus;          % GBPC

for i = 1:length(d_rand_sample)
      
        % Set time span
        tspan = 0:deltaT:final_time;
        
        % Preallocate Space
        RasBprime = zeros(length(tspan), 1);
        GBPCprime = zeros(length(tspan), 1);
        MCORprime = zeros(length(tspan), 1);
        MHCKAprime = zeros(length(tspan), 1);
        
        % Put parameter values into vector
        params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
        k1minus, k2minus, k3minus, k4minus, k5minus];
        
        % Calculate Steady State Values
        SSsolns = steadyStateSolutions(R, params);
        GBGSS =  SSsolns(:,1);  % Steady State Value of GBG
        
        % Use PDE for GBG
        GBG = pdefxn(final_time, d_rand_sample(i), lengthScale, deltaT, GBGSS, k1minus);
        
        % Set initial conditions
        init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
           SSsolns(:,5)];
        GBPCprime(1) = init_cond(1);
        MCORprime(1) = init_cond(2);
        RasBprime(1) = init_cond(3);
        MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
        % Calculate the percent change of myosin
        figure(4)
        plot(tspan, MCORprime, 'color', [.7 .7 .7], 'linewidth', 2)
        hold on
end


% Fixed Parameters
k2plus = init_alpha*k1plus;           % GBPC
k3plus = init_alpha*k1plus;           % RasB
k4plus = init_alpha*k1plus;           % MHCKA
k5plus = init_alpha*k1plus;           % MHCKA
k3minus = init_alpha*k1plus; 
     
% Set time span
tspan = 0:deltaT:final_time;

% Preallocate Space
RasBprime = zeros(length(tspan), 1);
GBPCprime = zeros(length(tspan), 1);
MCORprime = zeros(length(tspan), 1);
MHCKAprime = zeros(length(tspan), 1);

% Put parameter values into vector
params = [k1plus, k2plus, k3plus, k4plus, k5plus, ...
k1minus, k2minus, k3minus, k4minus, k5minus];

% Calculate Steady State Values
SSsolns = steadyStateSolutions(R, params);
GBGSS =  SSsolns(:,1);  % Steady State Value of GBG

% Use PDE for GBG
GBG = pdefxn(final_time, d, lengthScale, deltaT, GBGSS, k1minus);

% Set initial conditions
init_cond = [SSsolns(:,2) SSsolns(:,3) SSsolns(:,4)...
   SSsolns(:,5)];
GBPCprime(1) = init_cond(1);
MCORprime(1) = init_cond(2);
RasBprime(1) = init_cond(3);
MHCKAprime(1) = init_cond(4);
        
        for jj = 1:length(tspan)-1    
            GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) - ...
                k2minus*GBPCprime(jj) - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj) - k6minus*MCORprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                k4minus*RasBprime(jj) - k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                 - k5minus*MHCKAprime(jj) - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
% Calculate the percent change of myosin
plot(tspan, MCORprime, 'color', 'm', 'linewidth', 2);
xlabel('\bf Time (Seconds)', 'FontSize', 17)
ylabel('\bf Concentration', 'FontSize', 17)

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/SAshadedPlots', figureName, '.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder



