function Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha, MACH, Clalpha)
%%
% Outputs:
%   Cd0       -> Coeficiente de arrasto m�nimo (sem interfer�ncia do AoA)
% Inputs:
%   Cd_corpo  -> Coeficiente de arrasto do corpo
%   Cd_base   -> Coeficiente de arrasto base
%   Cd_aleta  -> Coeficiente de arrasto das aletas
%   Cd_interf -> Coeficiente de arrasto de interfer�ncia
%   Va        -> M�dulo da velocidade
%   Geral     -> Estrutura com dados gerais do minifoguete
%   Aleta     -> Estrutura com dados das aletas do minifoguete
%   alpha     -> �ngulo de ataque em radianos
%   MACH      -> N�mero de MACH
%   Clalpha   -> Coeficiente de sustenta��o
%%

% Reynolds de transiss�o cr�tico
ReTC = 5e5;
% N�mero de Reynolds do corpo
Re_corpo = Va * Geral.Ltr / Geral.ViscAr;
% N�mero de Reynolds das aletas
Re_aleta = Va * Aleta.Lm / Geral.ViscAr;
% Par�metro B de Reynolds
B_corpo = ReTC * ((0.074/Re_corpo^(1/5)) - (1.328/sqrt(Re_corpo)));
B_aleta = ReTC * ((0.074/Re_aleta^(1/5)) - (1.328/sqrt(Re_aleta)));

if Re_corpo < ReTC
    Cf_corpo = 1.328 / sqrt(Re_corpo);
else
    Cf_corpo = (0.074/Re_corpo^(1/5)) - (B_corpo/Re_corpo);
end

if Re_aleta < ReTC
    Cf_aleta = 1.328 / sqrt(Re_aleta);
else
    Cf_aleta = (0.074/Re_aleta^(1/5)) - (B_aleta/Re_aleta);
end

% Coeficiente de arrasto m�nimo
Cd0 = Cd_corpo * Cf_corpo + Cd_base + Cd_aleta * Cf_aleta + Cd_interf * Cf_aleta;

% Ajuste com o coeficiente de sustenta��o
% Cd0 = (Cd0*cos(alpha) - 0.5*Clalpha*sin(2*alpha)) / (1 - sin(alpha)^2);

% Ajuste para fluido incompresss�vel
Cd0 = Cd0 / sqrt(1 - MACH^2);
