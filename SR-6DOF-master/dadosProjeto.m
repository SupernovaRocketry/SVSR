function [Prop, Geral, Coifa, Corpo, Aleta, Paraq] = dadosProjeto()
%%
% Outputs:
%   Prop  -> Estrutura com os dados do motor foguete
%   Geral -> Estrutura com os dados gerais do minifoguete
%   Coifa -> Estrutura com os dados da coifa do minifoguete
%   Corpo -> Estrutura com os dados do corpo do minifoguete
%   Aleta -> Estrutura com os dados das aletas do minifoguete
%%

% Dados de Propulsão
Prop.Thrust_N = xlsread('Calculos(Nakka).xls','Performance','J28:J910');    % Thrust (N)
Prop.Thrust_t = xlsread('Calculos(Nakka).xls','Performance','L28:L910');    % Tempo (s)
Prop.ImpulsoTotal_N = 729;                                                  % Impulso total (N.s)
Prop.MassaPropelente = 0.569;                                               % Massa do propelente (kg)
Prop.MassaTotalMotor = 0.820 + Prop.MassaPropelente;                        % Massa do motor (kg)

% Dados Gerais
Geral.Ltr = 1.67757;                                                    % Comprimento total do foguete (m)
Geral.m = 1500e-3;                                                      % Massa do foguete vazio [corpo + aletas + coifa] (kg)
Geral.MassaEletronica = 800e-3;                                         % Massa da eletrônica (kg)
Geral.mtotal = Geral.m + Prop.MassaTotalMotor + Geral.MassaEletronica;  % Massa total do foguete [Massa do foguete vazio + Massa do motor foguete + Massa da eletrônica] (kg)
Geral.ViscAr = 15.11e-6;                                                % Viscosidade cinemática do ar seco (m^2/s) - considerando uma temperatura constante de 20°C
Geral.xcg = 1080e-3;                                                    % Posição do CG medida a partir da ponta da coifa (m)
Geral.xcp = 1220e-3;                                                    % Posição do CP medida a partir da ponta da coifa (m)

% Dados do Corpo
Corpo.Lb = 1.3;                             % Comprimento do corpo (m)
Corpo.Lc = 40e-3;                           % Comprimento da saia (m)
Corpo.Db = 67.5e-3;                         % Diâmetro do foguete (m)
Corpo.Du = Corpo.Db;                        % Diâmetro no começo da saia (m)
Corpo.Dd = 45e-3;                           % Diâmetro no final da saia
Corpo.AreaFoguete = (Corpo.Db)^2*pi/4;      % Seção do foguete (m^2)
Corpo.EspessuraFoguete = 2e-3;              % Espessura do corpo (m)
Corpo.SurfFoguete = 275674.76e-6;           % Superfície de contato do corpo (m^2)
Corpo.RugosidadeFoguete = 8;                % Rugosidade do corpo (micrometro)
Corpo.K = 0.5;                              % Constante definida no projeto do corpo

% Dados das Aletas
Aleta.n = 3;                                                % Número de aletas
Aleta.Lr = 160e-3;                                          % Comprimento da base da aleta (m)
Aleta.Lt = 28.82e-3;                                        % Comprimento da ponta da aleta (m)
Aleta.Lm = 110.11e-3;                                       % Comprimento médio da aleta (m)
Aleta.Ls = 72.8e-3;                                         % Semi span da aleta(m)
Aleta.Lts = 179.2e-3;                                       % Comprimento do span total
Aleta.Df = Corpo.Db;                                        % Diâmetro do corpo no começo da aleta (m)
Aleta.EspessuraAleta = 1e-3;                                % Espessura das aletas (m)
Aleta.SurfAleta = 6567.24e-6;                               % Superfície de contato das aletas em uma face (m^2) - Medida Solidworks
Aleta.CoefConstr = 1 + 2 * (Aleta.EspessuraAleta/Aleta.Lm); % Coeficiente construtivo das aletas
Aleta.RugosidadeAleta = 5;                                  % Rugosidade da aleta (micrometro)
Aleta.Kf = 1 + ((Aleta.Df/2) / (Aleta.Ls + (Aleta.Df/2)));  % Constante definida no projeto da aleta

% Dados da Coifa
Coifa.Ln = 337.57e-3;           % Comprimento da coifa (m)
Coifa.SurfCoifa = 49615.96e-6;  % Superfície de contato da coifa (m^2) - Medida Solidworks
Coifa.Dn = Corpo.Db;            % Diâmetro máximo da Coifa (m)
Coifa.RugosidadeCoifa = 50;     % Rugosidade da coifa (micrometro)

% Dados dos Paraquedas
Paraq.R1 = 110e-3;              % Raio do paraquedas inicial (m)
Paraq.R2 = 200e-3;              % Raio do paraquedas principal (m)
Paraq.S1 = pi * Paraq.R1^2;     % Área da seção transversal do paraquedas inicial (m^2)
Paraq.S2 = pi * Paraq.R2^2;     % Área da seção transversal do paraquedas principal (m^2)
Paraq.Coef = 2;                 % Coeficiente aerodinâmico dos paraquedas
