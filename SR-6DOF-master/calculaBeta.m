function [deg,rad] = calculaBeta(in1,in2)
%%
% Inputs:
%   in1 -> Componente do numerador
%   in2 -> Componente do denominador
% Outputs:
%   deg -> �ngulo em graus (arco seno)
%   rad -> �ngulo em radianos (arco seno)
%%
% �ngulo em radianos
rad = asin(in1/in2);
% �ngulo em graus
deg = rad2deg(rad);
