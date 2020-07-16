function Cmalpha = calculaCoefMomento(Geral, Clalpha, cref)
%%
% Inputs:
%   Geral   -> Estrutura com os dados gerais do minifoguete
%   Clalpha -> Coeficiente de Sustenta��o
%   cef     -> Corda de refer�ncia (comprimento da base da aleta)
% Outputs:
%   Cmalpha -> Coeficiente de momentos aerodin�micos
%%
Cmalpha = (Geral.xcg - Geral.xcp) * Clalpha / cref;
