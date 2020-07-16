function [rho, P, Temp_F, Temp_C, g] = calculaDadosAtm(h)
%% 
% Inputs:
%   h       -> altitude base do referencial inercial em metros
% Outputs:
%   rho     -> Densidade do ar em rela��o ao referencial inercial
%   P       -> Press�o atmosf�rica em rela��o ao referencial inercial
%   Temp_F  -> Temperatura atmosf�rica em �C em rela��o ao ref. inercial
%   Temp_C  -> Temperatura atmosf�rica em �F
%   g       -> Acelera��o da gravidade
%% 

% C�lculo da densidade do ar - linear no intervalo de 0 a 2000m (rho = (-1.05e-4 * h) + 1.22)
rho_mar = 1.22;                 % Densidade do ar no n�vel do mar (kg/m^3)
rho = (-1.05e-4 * h) + rho_mar;

% Dados para c�lculo da press�o atmosf�rica pela altitude
P0 = 101.325;                           % Press�o atmosf�rica padr�o ao n�vel do mar (kPa)
g = 9.81;                               % Acelera��o da gravidade (m/s^2)
M = 0.02896968;                         % Massa molar do ar seco (kg/mol)
T0 = ((288.15 - 273.15)*(9/5)) + 32;    % Temperatura ao n�vel do mar (�F)
R0 = 8.31582991;                        % Constante do g�s universal (J/(mol.K))
P = P0 * exp(-(g * h * M) / (T0 * R0)); % Press�o atmosf�rica

% Dados para c�lculo da temperatura pela altitude
Tbase = 20;                         % Temperatura m�dia no ano na cidade de Juiz de Fora (�C)
LapseRate = 9.8 * 10^(-3);          % Lapse Rate - o quanto a temperatura cai com a altitude (�C/m)
Temp_C = -LapseRate * h + Tbase;    % Temperatura atmosf�rica em �C
Temp_F = (Temp_C * 9/5) + 32;       % Temperatura atmosf�rica em �F
