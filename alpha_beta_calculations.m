% Clear the workspace
clc
clear
close all

% Set paramets values
[R, lengthScale, d, deltaT, final_time,...
    k1plus, alpha, k2plus, k3plus, k4plus, k5plus, k3minus,...
    k1minus, beta, k2minus, k4minus, k5minus, k6minus,...
    k_0] = setParameters();

gamma = (k1minus/k1plus)*(1 + (k1minus/k1plus)*(1/R));
sigma = k1plus/k1minus;

max_beta = ((-1 + sqrt(1 + 4*gamma*sigma)))/(2*gamma)

