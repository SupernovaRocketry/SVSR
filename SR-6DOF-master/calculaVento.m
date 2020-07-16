function gamma = calculaVento(K, dir)
%%
% Outputs:
%   gamma -> Vetor de vento lateral no referencial inercial
% Inputs:
%   K     -> Intensidade do vento
%   dir   -> String com a direção do vento
%%
switch dir
    case 'norte'
        gamma = K * [-1 0 0 0 0 0]';
    case 'leste'
        gamma = K * [0 -1 0 0 0 0]';
    case 'sul'
        gamma = K * [1 0 0 0 0 0]';
    case 'oeste'
        gamma = K * [0 1 0 0 0 0]';
    case 'nordeste'
        gamma = K * [-0.5 -0.5 0 0 0 0]';
    case 'sudeste'
        gamma = K * [0.5 -0.5 0 0 0 0]';
    case 'sudoeste'
        gamma = K * [0.5 0.5 0 0 0 0]';
    case 'noroeste'
        gamma = K * [-0.5 0.5 0 0 0 0]';
end
