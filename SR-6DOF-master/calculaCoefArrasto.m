function [Cd_corpo, Cd_aleta, Cd_interf, Cd_base] = calculaCoefArrasto(Geral, Coifa, Corpo, Aleta)
%%
% Inputs:
%   Geral     -> Estrutura com dados gerais do minifoguete
%   Coifa     -> Estrutura com dados da coifa do minifoguete
%   Corpo     -> Estrutura com dados do corpo do minifoguete
%   Aleta     -> Estrutura com dados das aletas do minifoguete
% Outputs:
%   Cd_corpo  -> Coeficiente de arrasto referente ao corpo
%   Cd_aleta  -> Cofeiciente de arrasto referente às aletas
%   Cd_interf -> Coeficiente de arrasto de interferência das aletas
%   Cd_base   -> Coeficiente de arrasto de base
%%

% Corpo
Cd_corpo = (1 + 60/(Geral.Ltr/Coifa.Dn)^3 + 0.0025*(Corpo.Lb/Corpo.Db)) * (2.7*(Coifa.Ln/Corpo.Db) + 4*(Corpo.Lb/Corpo.Db) + 2*(1-(Corpo.Dd/Corpo.Db))*(Corpo.Lc/Corpo.Db)); % Ainda muliplica o Cf do corpo

% Base
Cd_base = 0.029 * (Corpo.Dd/Coifa.Dn)^3 / sqrt(Cd_corpo);

% Aleta
Afe = 0.5 * (Aleta.Lr + Aleta.Lt) * Aleta.Ls;
Afp = Afe + 0.5 * Aleta.Df * Aleta.Lr;
Cd_aleta = 2 * Aleta.CoefConstr * (4 * Aleta.n * Afp)/(pi * Aleta.Df^2); % Ainda multiplica o Cf das aletas

% Interferência
Cd_interf = 2 * Aleta.CoefConstr * (4 * Aleta.n * (Afp - Afe)) / (pi * Aleta.Df^2); % Ainda multiplica o Cf das aletas
