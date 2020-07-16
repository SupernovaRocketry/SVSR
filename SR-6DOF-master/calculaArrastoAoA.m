function Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad)
%%
% Output: 
%   Cdi       -> Coeficiente de arrasto de interferência (devido ao AoA)
% Input:
%   Geral     -> Estrutura com dados gerais do minifoguete
%   Coifa     -> Estrutura com dados da coifa do minifoguete
%   Corpo     -> Estrutura com dados do corpo do minifoguete
%   Aleta     -> Estrutura com dados das aletas do minifoguete
%   alpha_rad -> Ângulo de ataque em radianos
%%

% Ajuste do AoA
Rs = Aleta.Lts / Aleta.Df;
Kfb = 0.8065*Rs^2 + 1.1553*Rs;
Kbf = 0.1935*Rs^2 + 0.8174*Rs + 1;
Afe = 0.5 * (Aleta.Lr + Aleta.Lt) * Aleta.Ls;
Afp = Afe + 0.5 * Aleta.Df * Aleta.Lr;
C1 = 1.2 * ((4 * Afp) / (pi * Aleta.Df^2));
C2 = 3.12 * (Kfb + Kbf -1) * ((4 * Afe) / (pi * Aleta.Df^2));

if (alpha_rad ~= 0)
    % Parâmetros de alpha
    delta = alpha_rad / sqrt(11.65 + alpha_rad^2);
    eta = 0.4482 + 0.1041*log(alpha_rad);

    % Ajuste do corpo
    Cdf = 2*delta*alpha_rad^2 + (3.6*eta*(1.36*Geral.Ltr - 0.55*Coifa.Ln)/(pi*Corpo.Db)) * alpha_rad^3;

    % Ajuste das aletas
    Cda = alpha_rad^2 * (C1 + C2);

    % Coeficiente de arrasto de interferência
    Cdi = Cdf + Cda;
else
    Cdi = 0;
end
