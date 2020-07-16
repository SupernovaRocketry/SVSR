function [Cl_coifa, Cl_corpo, Cl_aleta, Cl_trans] = calculaCoefSustentacao(Geral, Coifa, Corpo, Aleta)
%%
% Inputs:
%   Geral    -> Estrutura com os dados gerais do minifoguete
%   Coifa    -> Estrutura com os dados da coifa do minifoguete
%   Corpo    -> Estrutura com os dados do corpo do minifoguete
%   Aleta    -> Estrutura com os dados das aletas do minifoguete
% Outputs:
%   Cl_coifa -> Coeficiente de sustentação da coifa
%   Cl_trans -> Coeficiente de sustentação de transições no corpo
%   Cl_corpo -> Coeficiente de sustentação do corpo
%   Cl_aleta -> Coeficiente de sustentação das aletas
%%

% Coifa
Cl_coifa = 2;

% Transição
Cl_trans = 2 * ((Corpo.Dd / Coifa.Dn)^2 - (Corpo.Du / Coifa.Dn)^2);

% Corpo
Cl_corpo = Corpo.K * (Corpo.SurfFoguete / Corpo.AreaFoguete); % (ainda tem o AoA multiplicando)

% Aleta
Cl_aleta = Aleta.Kf * ((4 * Aleta.n * (Aleta.Ls / Coifa.Dn)^2) / (1 + sqrt(1 + (2*Aleta.Lm /(Aleta.Lr + Aleta.Lt)))));
