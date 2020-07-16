function dvdt = edo6DOF(t,v,M,C,D,G,Thrust_t,Thrust_N)
%% Fun��o original do modelo:
% 
% M * dvdt + C * v + D * v + G = tau
% 
%% Interpola��o plo tempo da EDO
Thrust = interp1(Thrust_t , Thrust_N , t);
%% Matriz de for�as e momentos
tau = [ Thrust
          0
          0
          0
          0
          0    ];
%% Representa��o da EDO em MATLAB
dvdt =  -(C/M) * v - (D/M) * v - (G/M) + tau;
