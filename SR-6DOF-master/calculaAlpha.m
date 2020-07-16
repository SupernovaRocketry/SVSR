function [deg,rad] = calculaAlpha(in1,in2)
%%
% Inputs:
%   in1 -> Componente do numerador
%   in2 -> Componente do denominador
% Outputs:
%   deg -> �ngulo em graus (arco tangente)
%   rad -> �ngulo em radianos (arco tangente)
%%
% �ngulo em radianos
rad = atan(in1/in2);
% �ngulo em graus
deg = rad2deg(rad);
