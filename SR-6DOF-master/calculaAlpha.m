function [deg,rad] = calculaAlpha(in1,in2)
%%
% Inputs:
%   in1 -> Componente do numerador
%   in2 -> Componente do denominador
% Outputs:
%   deg -> Ângulo em graus (arco tangente)
%   rad -> Ângulo em radianos (arco tangente)
%%
% Ângulo em radianos
rad = atan(in1/in2);
% Ângulo em graus
deg = rad2deg(rad);
