function [out] = calculaModulo(in1,in2,in3)
%%
% Inputs:
%   in1 -> Componente em x
%   in2 -> Componente em y
%   in3 -> Componente em z
% Outputs:
%   out -> Módulo das 3 componentes
%%
out = sqrt(in1^2 + in2^2 + in3^2);
