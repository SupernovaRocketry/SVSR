close all
clear all
clc
%% =======================================================================
inicio = tic;
%% Dados de Projeto

[Prop, Geral, Coifa, Corpo, Aleta, Paraq] = dadosProjeto();

%% Coeficientes de Arrasto

[Cd_corpo, Cd_aleta, Cd_interf, Cd_base] = calculaCoefArrasto(Geral, Coifa, Corpo, Aleta);

%% Coeficientes de Sustentação

[Cl_coifa, Cl_corpo, Cl_aleta, Cl_trans] = calculaCoefSustentacao(Geral, Coifa, Corpo, Aleta);

%% Dados Atmosféricos

% Altitude base do referencial inercial em metros
h = 800;

[rho, P, Temp_F, Temp_C, g] = calculaDadosAtm(h);

%% *******************************************************
%% Modelo Estocástico
%% *******************************************************

% Número de simulações
n = 2;

for k = 1:n

clear eta eta_i eta_int eta_p1 eta_p2 eta_ponto eta_ponto_p1 eta_ponto_p2 eta_pontoX eta_pontoY eta_pontoZ etaX etaY etaZ az_p1 az_p2 v t t_p1 t_p2

Geral.mtot = Geral.mtotal + ((0.3)*randn);

%% Momentos de Inércia

Ix = 687989740.24e-9 + ((0.1)*randn);
Iy = 812188074.87e-9 + ((0.1)*randn);
Iz = 126562111.96e-9 + ((0.1)*randn);
Ixy = 667702.56e-9;
Ixz = 49300684.86e-9;
Iyz = 224734.27e-9;
Iyx = Ixy;
Izx = Ixz;
Izy = Iyz;

I = [ Ix   -Ixy  -Ixz
     -Iyx   Iy   -Iyz
     -Izx  -Izy   Iz  ];

%% *******************************************************
%% Inicialização do Simulador
%% *******************************************************

% Matriz de posição (Referencial inercial)
x = 0;                      % 
y = 0;                      % POSIÇÕES DE LANÇAMENTO
z = 0;                      % 
Langle = 85 + (89-85)*rand; % Ângulo de lançamento
phi_x = 0;                  % Roll ou rolamento
phi_y = deg2rad(Langle);    % Pitch ou arfagem
phi_z = 0;                  % Yaw ou guinada

eta_i = [ x
          y
          z
          phi_x
          phi_y
          phi_z ];

% Matriz de velocidades (Referencial do corpo)
vx = 1;     %
vy = 0;     % VELOCIDADES LINEARES
vz = 0;     %
wx = 0;     %
wy = 0;     % VELOCIDADES ANGULARES
wz = 0;     %

v_i = [ vx
        vy
        vz
        wx
        wy
        wz ];

% Módulo da velocidade
Va = calculaModulo(vx,vy,vz);
% Número de MACH
MACH = Va / 360;

% Matriz de forças (Referencial do corpo)
Fx = 0;     %
Fy = 0;     % FORÇAS
Fz = 0;     %
Mx = 0;     %
My = 0;     % MOMENTOS
Mz = 0;     %

tau = [ Fx
        Fy
        Fz
        Mx
        My
        Mz ];

% Ângulos de ataque e escorregamento iniciais
[alpha_deg, alpha_rad] = calculaAlpha(vz,vx);
[beta_deg, beta_rad] = calculaBeta(vy,Va);

% Coeficiente de Sustentação (Clalpha)
Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans + ((2)*randn);

% Coeficiente de Arrasto (Cd0)
Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha) + ((0.1)*randn);
Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad) + ((0.3)*randn);

% Coeficiente de Momento Aerodinâmico (Cmalpha)
cref = Aleta.Lr;
Cmalpha = calculaCoefMomento(Geral, Clalpha, cref) + ((2)*randn);

% Matriz de transformação de coordenadas (Corpo > Inercial)
J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
      sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
     -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
      0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
      0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
      0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

%% Dados de Vento Lateral

K = 1 + 8*rand;
possiveldir = ["norte","sul","leste","oeste","sudeste","sudoeste","nordeste","noroeste"];
gamma = calculaVento(K, possiveldir(randi(numel(possiveldir))));
gammar = J\gamma;

% Matriz de momentos
M = [ Geral.mtot  0           0            0     0    0
      0           Geral.mtot  0            0     0    0
      0           0           Geral.mtot   0     0    0
      0           0           0            Ix   -Ixy  -Ixz
      0           0           0           -Iyx   Iy   -Iyz
      0           0           0           -Izx  -Izy   Iz  ];

% Matriz dos termo de Coriolis
C = [ 0            0            0                   0              Geral.mtot*vz     -Geral.mtot*vy
      0            0            0             -Geral.mtot*vz            0             Geral.mtot*vx
      0            0            0              Geral.mtot*vy      -Geral.mtot*vx           0
      0            0            0                   0                   0                  0
      0            0            0                -Iy*wz                 0                Ix*wx
      0            0            0                 Iy*wy              -Ix*wx                0        ];

% Matriz de forças aerodinâmicas
D = 0.5 * rho * Corpo.AreaFoguete * [ Cd0*Va   (Cdi-Clalpha)*vy          (Cdi-Clalpha)*vz          0  0  0
                                      Cd0*vy   Clalpha*Va+(Cdi*vy^2/Va)  Cdi*vz*vy/Va              0  0  0
                                      Cd0*vz   Cdi*vy*vz/Va              Clalpha*Va+(Cdi*vz^2/Va)  0  0  0
                                      0        0                         0                         0  0  0
                                      0        0                        -Cmalpha*cref*Va           0  0  0
                                      0        Cmalpha*cref*Va           0                         0  0  0 ];

% Matriz de força gravitacional
G = -Geral.mtot * g * [    -sin(phi_y)
                       cos(phi_y)*sin(phi_x)
                       cos(phi_y)*cos(phi_x)
                                0
                                0
                                0            ];

% Velocidade inicial (referencial inercial)
eta_ponto(:,1) = J * v_i;

%% *******************************************************
%% Resolução do Método
%% *******************************************************

imax = 5000;         % Valor máximo da variável de controle
dt = 0.01;           % Passo da simulação
t(1) = 0;            % Condição inicial de tempo
v(:,1) = v_i;        % Condição inicial de velocidade

% Interpolação do expuxo
Thrust = interp1(Prop.Thrust_t , Prop.Thrust_N , 0:dt:Prop.Thrust_t(end));
Thrust = Thrust + (50*randn(1,length(Thrust)));
Thrust = [Thrust , zeros(1,1e6)];

% EDO do sistema
F_xy = @(tempo,vel,tau) -(M\C) * vel - (M\D) * vel - (M\G) + (M\tau);

%% Forças atuantes durante lançamento
for i = 1:30
    % Matrizes de forças
    D = zeros(6,6);
    tauSim = [ Thrust(i+1)
                  0
                  0
                  0
                  0
                  0        ];
    %% Método de Runge-Kutta
    % Coeficientes do equacionamento de Runge-Kutta
    k_1 = F_xy(t(i) , v(:,i) , tauSim);
    k_2 = F_xy(t(i)+0.5*dt , v(:,i)+0.5*dt*k_1 , tauSim);
    k_3 = F_xy((t(i)+0.5*dt) , (v(:,i)+0.5*dt*k_2) , tauSim);
    k_4 = F_xy((t(i)+dt) , (v(:,i)+k_3*dt) , tauSim);
    % Equação principal de Runge-Kutta
    v(:,i+1) = v(:,i) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*dt;
    
    %% Atualização das variáveis da simulação
    % Incremento do tempo
    t(i+1) = t(i) + dt;

    % Condição da base de lançamento (velocidades perpendiculares = 0 ; velocidades angulares = 0)
    v(2,i+1) = 0;
    v(3,i+1) = 0;
    v(4,i+1) = 0;
    v(5,i+1) = 0;
    v(6,i+1) = 0;

    % Vetor velocidade no triedro inercial
    eta_ponto(:,i+1) = J * v(:,i+1) - gammar;

    % Vetor posição no triedro inercial
    eta_int = cumtrapz(t , eta_ponto , 2);
    eta(:,i+1) = eta_int(:,i+1) + eta_i;

    % Simplificações de coordenadas
    vx = v(1,i+1);
    vy = v(2,i+1);
    vz = v(3,i+1);
    wx = v(4,i+1);
    wy = v(5,i+1);
    wz = v(6,i+1);

    phi_x = eta(4,i+1);
    phi_y = deg2rad(Langle);
    phi_z = eta(6,i+1);

    % Módulo da velocidade
    Va = calculaModulo(vx,vy,vz);
    % Número de MACH
    MACH = Va / 360;

    % Ângulos de ataque e escorregamento
    [alpha_deg, alpha_rad] = calculaAlpha(vz,vx);
    [beta_deg, beta_rad] = calculaBeta(vy,Va);

    % Coeficiente de Sustentação (Clalpha)
    Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans;

    % Coeficiente de Arrasto (Cd0)
    Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha);
    Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad);

    % Coeficiente de Momento Aerodinâmico (Cmalpha)
    cref = Aleta.Lr;
    Cmalpha = calculaCoefMomento(Geral, Clalpha, cref);

    %% Atualização das matrizes
    % Matriz de transformação de coordenadas (Corpo > Inercial)
    J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
          sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
         -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
          0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
          0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
          0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

    % Matriz de momentos
    M = [ Geral.mtot  0           0            0     0    0
          0           Geral.mtot  0            0     0    0
          0           0           Geral.mtot   0     0    0
          0           0           0            Ix   -Ixy  -Ixz
          0           0           0           -Iyx   Iy   -Iyz
          0           0           0           -Izx  -Izy   Iz  ];

    % Matriz de Coriolis
    C = [ 0            0            0                   0             Geral.mtot*vz     -Geral.mtot*vy
          0            0            0            -Geral.mtot*vz             0            Geral.mtot*vx
          0            0            0             Geral.mtot*vy      -Geral.mtot*vx            0
          0            0            0                   0                   0                  0
          0            0            0                -Iy*wz                 0                Ix*wx
          0            0            0                 Iy*wy              -Ix*wx                0        ];

    % Matriz de força gravitacional
    G = -Geral.mtot * g * [    -sin(phi_y)
                           cos(phi_y)*sin(phi_x)
                           cos(phi_y)*cos(phi_x)
                                    0
                                    0
                                    0            ];
end
i = i + 1;

%% Forças atuantes no espaço livre
while true
    % Definição da matriz tau
    tauSim = [ Thrust(i+1)
                  0
                  0
                  0
                  0
                  0        ];
    %% Método de Runge-Kutta
    if (eta_ponto(3,i) > 0.1 || i == imax)  % Se vz > 0.1 (descida) ou i = imax (overflow)
        break;
    else
        % Coeficientes do equacionamento de Runge-Kutta
        k_1 = F_xy(t(i) , v(:,i) , tauSim);
        k_2 = F_xy(t(i)+0.5*dt , v(:,i)+0.5*dt*k_1 , tauSim);
        k_3 = F_xy((t(i)+0.5*dt) , (v(:,i)+0.5*dt*k_2) , tauSim);
        k_4 = F_xy((t(i)+dt) , (v(:,i)+k_3*dt) , tauSim);
        % Equação principal de Runge-Kutta
        v(:,i+1) = v(:,i) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*dt;
    end
    %% Atualização das variáveis da simulação
    % Incremento do tempo
    t(i+1) = t(i) + dt;

    % Vetor velocidade no triedro inercial
    eta_ponto(:,i+1) = J * v(:,i+1) - gammar;

    % Vetor posição no triedro inercial
    eta_int = cumtrapz(t, eta_ponto , 2);
    eta(:,i+1) = eta_int(:,i+1) + eta_i;

    % Simplificações de coordenadas
    vx = v(1,i+1);
    vy = v(2,i+1);
    vz = v(3,i+1);
    wx = v(4,i+1);
    wy = v(5,i+1);
    wz = v(6,i+1);

    phi_x = eta(4,i+1);
    phi_y = eta(5,i+1);
    phi_z = eta(6,i+1);

    % Módulo da velocidade
    Va = calculaModulo(vx,vy,vz);
    % Número de MACH
    MACH = Va / 360;

    % Ângulos de ataque e escorregamento
    [alpha_deg, alpha_rad] = calculaAlpha(vz,vx);
    [beta_deg, beta_rad] = calculaBeta(vy,Va);

    % Coeficiente de Sustentação (Clalpha)
    Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans + ((2)*randn);

    % Coeficiente de Arrasto (Cd0)
    Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha) + ((0.1)*randn);
    Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad) + ((0.3)*randn);

    % Coeficiente de Momento Aerodinâmico (Cmalpha)
    cref = Aleta.Lr;
    Cmalpha = calculaCoefMomento(Geral, Clalpha, cref) + ((2)*randn);

    %% Atualização das matrizes
    % Matriz de transformação de coordenadas (Corpo > Inercial)
    J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
          sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
         -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
          0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
          0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
          0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

    % Matriz de momentos
    M = [ Geral.mtot  0           0            0     0    0
          0           Geral.mtot  0            0     0    0
          0           0           Geral.mtot   0     0    0
          0           0           0            Ix   -Ixy  -Ixz
          0           0           0           -Iyx   Iy   -Iyz
          0           0           0           -Izx  -Izy   Iz  ];

    % Matriz de Coriolis
    C = [ 0            0            0                   0              Geral.mtot*vz    -Geral.mtot*vy
          0            0            0            -Geral.mtot*vz              0           Geral.mtot*vx
          0            0            0             Geral.mtot*vy       -Geral.mtot*vx           0
          0            0            0                   0                   0                  0
          0            0            0                -Iy*wz                 0                Ix*wx
          0            0            0                 Iy*wy              -Ix*wx                0        ];

    % Matriz de forças aerodinâmicas
    D = 0.5 * rho * Corpo.AreaFoguete * [ Cd0*Va    (Cdi-Clalpha)*vy           (Cdi-Clalpha)*vz          0  0  0
                                          Cd0*vy    Clalpha*Va+(Cdi*vy^2/Va)   Cdi*vz*vy/Va              0  0  0
                                          Cd0*vz    Cdi*vy*vz/Va               Clalpha*Va+(Cdi*vz^2/Va)  0  0  0
                                          0         0                          0                         0  0  0
                                          0         0                         -Cmalpha*cref*Va           0  0  0
                                          0         Cmalpha*cref*Va            0                         0  0  0 ];

    % Matriz de força gravitacional
    G = -Geral.mtot * g * [    -sin(phi_y)
                           cos(phi_y)*sin(phi_x)
                           cos(phi_y)*cos(phi_x)
                                    0
                                    0
                                    0            ];

    %% Incremento da variável de controle
    i = i + 1;
end

%% Acionamento do paraquedas inicial

ii = 1;         % Variável de controle do paraquedas inicial
iimax = 10000;  % Valor máximo da variável de controle
t_p1(1) = 0;    % Condição inicial de tempo do paraquedas inicial

% Condições iniciais de velocidade (referencial inercial)
if (gamma(1) == 0)
    eta_pontoX = eta_ponto(1,end);
else
    eta_pontoX = gamma(1);
end
if (gamma(2) == 0)
    eta_pontoY = eta_ponto(2,end);
else
    eta_pontoY = gamma(2);
end
eta_pontoZ = eta_ponto(3,end);
eta_ponto_p1 = [ eta_pontoX
                 eta_pontoY
                 eta_pontoZ
                     0
                     0
                     0      ];

% Condições iniciais de posição (referencial inercial)
etaX = eta(1,end);
etaY = eta(2,end);
etaZ = eta(3,end);
eta_p1 = [ etaX
           etaY
           etaZ
            0
            0
            0   ];

while true
    if(eta_p1(3,ii) > -500 || ii == iimax)
        break;
    else
        % Incremento do tempo
%         t(i+1) = t(i) + dt;
        t_p1(ii+1) = t_p1(ii) + dt;
        
        %% Regime com o paraquedas acionado
        D = 0.5 * rho * Paraq.Coef * Paraq.S1 * eta_ponto_p1(3,ii)^2;
        az_p1(ii+1) = g - (D / (Geral.mtotal - Prop.MassaPropelente));
        
        % Vetor de velocidades (referencial inercial)
        eta_pontoX(ii+1) = eta_pontoX(ii);
        eta_pontoY(ii+1) = eta_pontoY(ii);
        eta_pontoZ = cumtrapz(t_p1, az_p1);
        
        eta_ponto_p1 = [ eta_pontoX
                         eta_pontoY
                         eta_pontoZ
                         zeros(1,ii+1)
                         zeros(1,ii+1)
                         zeros(1,ii+1) ];
        
        % Vetor de posições
        eta_p1 = [ etaX
                   etaY
                   etaZ
                    0
                    0
                    0    ] + cumtrapz(t_p1, eta_ponto_p1, 2);

        %% Incremento das variáveis de controle
%         i = i + 1;
        ii = ii + 1;
    end
end

%% Acionamento do paraquedas principal

iii = 1;         % Variável de controle do paraquedas inicial
iiimax = 10000;  % Valor máximo da variável de controle
t_p2(1) = 0;     % Condição inicial de tempo do paraquedas inicial

% Condições iniciais de velocidade (referencial inercial)
if (gamma(1) == 0)
    eta_pontoX = eta_ponto_p1(1,end);
else
    eta_pontoX = gamma(1);
end
if (gamma(2) == 0)
    eta_pontoY = eta_ponto_p1(2,end);
else
    eta_pontoY = gamma(2);
end
eta_pontoZ = eta_ponto_p1(3,end);
eta_ponto_p2 = [ eta_pontoX
                 eta_pontoY
                 eta_pontoZ
                     0
                     0
                     0      ];

% Condições iniciais de posição (referencial inercial)
etaX = eta_p1(1,end);
etaY = eta_p1(2,end);
etaZ = eta_p1(3,end);
eta_p2 = [ etaX
           etaY
           etaZ
            0
            0
            0   ];

while true
    if(eta_p2(3,iii) > -0.1 || iii == iiimax)
        break;
    else
        % Incremento do tempo
%         t(i+1) = t(i) + dt;
        t_p2(iii+1) = t_p2(iii) + dt;
        
        %% Regime com o paraquedas acionado
        D = 0.5 * rho * Paraq.Coef * (Paraq.S1 + Paraq.S2) * eta_ponto_p2(3,iii)^2;
        az_p2(iii+1) = g - (D / (Geral.mtotal - Prop.MassaPropelente));
        
        % Vetor de velocidades (referencial inercial)
        eta_pontoX(iii+1) = eta_pontoX(iii);
        eta_pontoY(iii+1) = eta_pontoY(iii);
        eta_pontoZ = cumtrapz(t_p2, az_p2);

        eta_ponto_p2 = [ eta_ponto_p1(1,end)    
                         eta_ponto_p1(2,end)
                         eta_ponto_p1(3,end)
                                 0
                                 0
                                 0           ] +  [ eta_pontoX
                                                    eta_pontoY
                                                    eta_pontoZ
                                                    zeros(1,iii+1)
                                                    zeros(1,iii+1)
                                                    zeros(1,iii+1) ];

        % Vetor de posições
        eta_p2 = [ etaX
                   etaY
                   etaZ
                    0
                    0
                    0   ] + cumtrapz(t_p2, eta_ponto_p2, 2);

        %% Incremento das variáveis de controle
%         i = i + 1;
        iii = iii + 1;
    end
end

% Outputs
outPosicao{k,1} = cat(2 , eta , eta_p1 , eta_p2);
outVelocidade{k,1} = cat(2 , eta_ponto , eta_ponto_p1 , eta_ponto_p2);
outLaunchAng(k,1) = Langle;
outGamma{k,1} = gamma;

end
toc(inicio)
return
%% Plots

% Deslocamento 3D
figure
set(gcf,'color','w')
% set(gcf,'Units','centimeters','Position',[8 8 25 14])
hold on
for i = 1:n
    plot3(outPosicao{i}(1,:) , outPosicao{i}(2,:) , -outPosicao{i}(3,:) , 'color' , [0 0.447 0.741] , 'Linewidth' , 1.5)
    grid on
end
zlim([0 2500])
box on
zlabel('z (m)')
ylabel('y (m)')
xlabel('x (m)')
view(-20, 30)

% Histograma dos apogeus
for i = 1:n
    hmax(i) = max(abs(outPosicao{i}(3,:)));
end
pdfh = normpdf(500:3000,mean(hmax),std(hmax));
figure
set(gcf,'color','w')
% set(gcf,'Units','centimeters','Position',[8 8 25 14])
histogram(hmax , 80 , 'normalization' , 'pdf')
hold on
plot(500:3000 , pdfh , 'Linewidth' , 2);
ylabel('PDF_{h_{max}}')
xlabel('h_{max}')
grid on
aa = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
bb = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend([aa bb] , ['\mu = ' num2str(mean(hmax))] , ['\sigma = ' num2str(std(hmax))])

% Histograma da região de aterrissagem
for i = 1:n
    Xfinal(i) = outPosicao{i}(1,end);
    Yfinal(i) = outPosicao{i}(2,end);
end
figure
set(gcf,'color','w')
% set(gcf,'Units','centimeters','Position',[8 8 25 14])
histogram2(Xfinal , Yfinal , 80 , 'DisplayStyle' , 'tile' , 'normalization' , 'pdf')
xlabel('x (m)')
ylabel('y (m)')
view(-38, 58)

% Histograma da distância euclidiana da distância de aterrissagem em
% relação à origem
dist = sqrt(Xfinal.^2 + Yfinal.^2);
figure
set(gcf,'color','w')
% set(gcf,'Units','centimeters','Position',[8 8 25 14])
histogram(dist , 80 , 'normalization' , 'pdf')
grid on
xlabel('Distância Euclidiana')
ylabel('PDF_{dist}')
