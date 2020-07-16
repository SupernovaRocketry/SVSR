function [deg,rad] = calculaBeta(in1,in2)
%%
% Inputs:
%   in1 -> Componente do numerador
%   in2 -> Componente do denominador
% Outputs:
%   deg -> Ângulo em graus (arco seno)
%   rad -> Ângulo em radianos (arco seno)
%%
% Ângulo em radianos
rad = asin(in1/in2);
% Ângulo em graus
deg = rad2deg(rad);
