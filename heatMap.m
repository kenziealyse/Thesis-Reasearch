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

% Set alpha and beta values for for loop
alpha_values = 0:1:20;   %alpha
beta_values = 0:1:20;    %beta

% Preallocate space
percentChange = zeros(1,length(beta_values));

% Set paramets values
R = 0.8;                         % Initial active receptors 
lengthScale = 2.725;             % Diffusion length scale
d = .011;                        % Diffusion Rate
deltaT = .01;                    % Time Step
final_time = 13;                 % Final Time
k1plus = 1/10;                   % GBgamma
k1minus = 1/120;                 % GBgamma  

for i = 1:length(alpha_values)
    alpha = alpha_values(i);                % alpha
    k2plus = alpha*k1plus;           % GBPC
    k3plus = alpha*k1plus;           % MCOR
    k4plus = alpha*k1plus;           % RasB
    k5plus = alpha*k1plus;           % MHCKA
    k3minus = alpha*k1plus;          % MCOR
    
    for p = 1:length(beta_values)
        beta = beta_values(p);                  % beta
        k2minus = beta*k1minus;          % GBPC
        k4minus = beta*k1minus;          % RasB
        k5minus = beta*k1minus;          % MHCKA
        k6minus = beta*k1minus;          % MCOR degradation
        
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
                 GBPCprime(jj+1) =  GBPCprime(jj) + (k2plus*GBGSS*(1 - GBPCprime(jj)) ...
                 - k3plus*GBPCprime(jj)*(1-MCORprime(jj)))*deltaT;  %eq (GBPC)
            MCORprime(jj+1) =  MCORprime(jj) + (k3plus*GBPCprime(jj)*(1 - MCORprime(jj)) -...
                k3minus*MCORprime(jj)*MHCKAprime(jj))*deltaT;  %eq (XMCOR)
            RasBprime(jj+1) = RasBprime(jj) + (k4plus*GBG(jj)*(1 - RasBprime(jj)) - ...
                 k5plus*RasBprime(jj)*(1-MHCKAprime(jj)))*deltaT; % eq (RASB)
            MHCKAprime(jj+1) = MHCKAprime(jj) + (k5plus*RasBprime(jj)*(1 - MHCKAprime(jj))...
                  - k3minus*MHCKAprime(jj)*MCORprime(jj))*deltaT;  %eq (XMCHKA)
        end
    
        % Calculate the percent change of myosin
        percentChange(p,i) = 100*((max(MCORprime) - MCORprime(1))/MCORprime(1));

    end 
end

figure(1)
imagesc(alpha_values, beta_values, percentChange);
xlabel('\bf \alpha Value','FontSize',17);
ylabel('\bf \beta Value','FontSize',17);
colorbar;
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/parameter','heatmap.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder


figure(2)
imagesc(alpha_values, beta_values, percentChange);
xlabel('\bf \alpha Value','FontSize',17);
ylabel('\bf \beta Value','FontSize',17);
colorbar;
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
hold on
plot(15,5, '.r', 'MarkerSize',20)

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

figure_name = ['/parameter','heatmapwithalphanadbeta.pdf'];   
dirPath = strcat('/','figures', figure_name); % Directory Path
saveas(gcf,[pwd dirPath]); % Save Figure in Folder