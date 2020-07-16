function [rho, P, Temp_F, Temp_C, g] = calculaDadosAtm(h)
%% 
% Inputs:
%   h       -> altitude base do referencial inercial em metros
% Outputs:
%   rho     -> Densidade do ar em relação ao referencial inercial
%   P       -> Pressão atmosférica em relação ao referencial inercial
%   Temp_F  -> Temperatura atmosférica em °C em relação ao ref. inercial
%   Temp_C  -> Temperatura atmosférica em °F
%   g       -> Aceleração da gravidade
%% 

% Cálculo da densidade do ar - linear no intervalo de 0 a 2000m (rho = (-1.05e-4 * h) + 1.22)
rho_mar = 1.22;                 % Densidade do ar no nível do mar (kg/m^3)
rho = (-1.05e-4 * h) + rho_mar;

% Dados para cálculo da pressão atmosférica pela altitude
P0 = 101.325;                           % Pressão atmosférica padrão ao nível do mar (kPa)
g = 9.81;                               % Aceleração da gravidade (m/s^2)
M = 0.02896968;                         % Massa molar do ar seco (kg/mol)
T0 = ((288.15 - 273.15)*(9/5)) + 32;    % Temperatura ao nível do mar (°F)
R0 = 8.31582991;                        % Constante do gás universal (J/(mol.K))
P = P0 * exp(-(g * h * M) / (T0 * R0)); % Pressão atmosférica

% Dados para cálculo da temperatura pela altitude
Tbase = 20;                         % Temperatura média no ano na cidade de Juiz de Fora (°C)
LapseRate = 9.8 * 10^(-3);          % Lapse Rate - o quanto a temperatura cai com a altitude (ºC/m)
Temp_C = -LapseRate * h + Tbase;    % Temperatura atmosférica em °C
Temp_F = (Temp_C * 9/5) + 32;       % Temperatura atmosférica em °F
