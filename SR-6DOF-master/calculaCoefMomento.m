function Cmalpha = calculaCoefMomento(Geral, Clalpha, cref)
%%
% Inputs:
%   Geral   -> Estrutura com os dados gerais do minifoguete
%   Clalpha -> Coeficiente de Sustentação
%   cef     -> Corda de referência (comprimento da base da aleta)
% Outputs:
%   Cmalpha -> Coeficiente de momentos aerodinâmicos
%%
Cmalpha = (Geral.xcg - Geral.xcp) * Clalpha / cref;
