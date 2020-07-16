close all
clear all
clc
%% =======================================================================
%% CÓDIGO DE TRAJETÓRIA COM 6DOF DA EQUIPE SUPERNOVA ROCKETRY
% 
%% Autores: Setor de Simulações e Pesquisa (Supernova Rocketry)
%           - Yan Furtado Coutinho
%           - Michel Bernardino Machado
%           - Amanda Rodrigues de Melo
%           - Daniel Moraes Barbosa
%% Data: XX/XX/20XX
%%
%   Esse código foi desenvolvido com o intuito de calcular com maior
% precisão a trajetória do minifoguete e todas as grandezes envolvidas
% no vôo com 6 graus de liberdade (6DOF).
%   Os dados do motor foguete são obtidos através da tabela SRM fornecida
% pelo setor de propulsão.
%   Os dados de coeficientes de sustentação e de arrasto e a loaclização
% do centro de pressão (CP) são calculados segundo o método de Barrowman.
%% =======================================================================
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

%% Dados de Vento Lateral

K = 10;
dir = 'norte';

gamma = calculaVento(K, dir);

%% Momentos de Inércia

Ix = 687989740.24e-9;
Iy = 812188074.87e-9;
Iz = 126562111.96e-9;
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
Langle = 88;                % Ângulo de lançamento
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
Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans;

% Coeficiente de Arrasto (Cd0)
Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha);
Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad);

% Coeficiente de Momento Aerodinâmico (Cmalpha)
cref = Aleta.Lr;
Cmalpha = calculaCoefMomento(Geral, Clalpha, cref);

% Matriz de transformação de coordenadas (Corpo > Inercial)
J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
      sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
     -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
      0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
      0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
      0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

% Matriz de momentos
M = [ Geral.mtotal  0             0              0     0    0
      0             Geral.mtotal  0              0     0    0
      0             0             Geral.mtotal   0     0    0
      0             0             0              Ix   -Ixy  -Ixz
      0             0             0             -Iyx   Iy   -Iyz
      0             0             0             -Izx  -Izy   Iz  ];

% Matriz dos termo de Coriolis
C = [ 0            0            0                    0               Geral.mtotal*vz     -Geral.mtotal*vy
      0            0            0             -Geral.mtotal*vz             0              Geral.mtotal*vx
      0            0            0              Geral.mtotal*vy      -Geral.mtotal*vx            0
      0            0            0                    0                     0                    0
      0            0            0                 -Iy*wz                   0                  Ix*wx
      0            0            0                  Iy*wy                -Ix*wx                  0          ];

% Matriz de forças aerodinâmicas
D = 0.5 * rho * Corpo.AreaFoguete * [ Cd0*Va   (Cdi-Clalpha)*vy          (Cdi-Clalpha)*vz          0  0  0
                                      Cd0*vy   Clalpha*Va+(Cdi*vy^2/Va)  Cdi*vz*vy/Va              0  0  0
                                      Cd0*vz   Cdi*vy*vz/Va              Clalpha*Va+(Cdi*vz^2/Va)  0  0  0
                                      0        0                         0                         0  0  0
                                      0        0                        -Cmalpha*cref*Va           0  0  0
                                      0        Cmalpha*cref*Va           0                         0  0  0 ];

% Matriz de força gravitacional
G = -Geral.mtotal * g * [    -sin(phi_y)
                         cos(phi_y)*sin(phi_x)
                         cos(phi_y)*cos(phi_x)
                                  0
                                  0
                                  0            ];

% Velocidade inicial (referencial inercial)
eta_ponto = J * v_i;

% Vento inicial (referencial do corpo)
gammar = J\gamma;

%% *******************************************************
%% Resolução do Método
%% *******************************************************

imax = 5000;         % Valor máximo da variável de controle
dt = 0.01;           % Passo da simulação
t(1) = 0;            % Condição inicial de tempo
v(:,1) = v_i;        % Condição inicial de velocidade

% Interpolação do expuxo
Thrust = [interp1(Prop.Thrust_t , Prop.Thrust_N , 0:dt:Prop.Thrust_t(end)) , zeros(1,1e6)];

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
    eta_ponto(:,i+1) = J * v(:,i+1);

    % Vetor posição no triedro inercial
    eta_int = cumtrapz(t, eta_ponto, 2);
    eta = eta_i + eta_int;

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
    MACH(i+1) = Va / 360;

    % Ângulos de ataque e escorregamento
    [alpha_deg(i+1), alpha_rad] = calculaAlpha(vz,vx);
    [beta_deg(i+1), beta_rad] = calculaBeta(vy,Va);

    % Coeficiente de Sustentação (Clalpha)
    Clalpha(i+1) = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans;

    % Coeficientes de Arrasto (Cd0 e Cdi)
    Cd0(i+1) = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH(i+1), Clalpha(i+1));
    Cdi(i+1) = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad);

    % Coeficiente de Momento Aerodinâmico (Cmalpha)
    cref = Aleta.Lr;
    Cmalpha(i+1) = calculaCoefMomento(Geral, Clalpha(i+1), cref);

    %% Atualização das matrizes
    % Matriz de transformação de coordenadas (Corpo > Inercial)
    J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
          sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
         -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
          0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
          0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
          0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

    % Matriz de momentos
    M = [ Geral.mtotal  0             0              0     0    0
          0             Geral.mtotal  0              0     0    0
          0             0             Geral.mtotal   0     0    0
          0             0             0              Ix   -Ixy  -Ixz
          0             0             0             -Iyx   Iy   -Iyz
          0             0             0             -Izx  -Izy   Iz  ];

    % Matriz de Coriolis
    C = [ 0            0            0                    0               Geral.mtotal*vz     -Geral.mtotal*vy
          0            0            0             -Geral.mtotal*vz             0              Geral.mtotal*vx
          0            0            0              Geral.mtotal*vy      -Geral.mtotal*vx            0
          0            0            0                    0                     0                    0
          0            0            0                 -Iy*wz                   0                  Ix*wx
          0            0            0                  Iy*wy                -Ix*wx                  0          ];

    % Matriz de força gravitacional
    G = -Geral.mtotal * g * [    -sin(phi_y)
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
    eta_ponto(:,i+1) = J * v(:,i+1);

    % Vetor posição no triedro inercial
    eta_int = cumtrapz(t, eta_ponto, 2);
    eta = eta_i + eta_int;

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
    MACH(i+1) = Va / 360;

    % Ângulos de ataque e escorregamento
    [alpha_deg(i+1), alpha_rad] = calculaAlpha(vz,vx);
    [beta_deg(i+1), beta_rad] = calculaBeta(vy,Va);

    % Coeficiente de Sustentação (Clalpha)
    Clalpha(i+1) = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans;

    % Coeficientes de Arrasto (Cd0 e Cdi)
    Cd0(i+1) = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH(i+1), Clalpha(i+1));
    Cdi(i+1) = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad);

    % Coeficiente de Momento Aerodinâmico (Cmalpha)
    cref = Aleta.Lr;
    Cmalpha(i+1) = calculaCoefMomento(Geral, Clalpha(i+1), cref);

    %% Atualização das matrizes
    % Matriz de transformação de coordenadas (Corpo > Inercial)
    J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
          sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
         -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
          0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
          0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
          0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

    % Matriz de momentos
    M = [ Geral.mtotal  0             0              0     0    0
          0             Geral.mtotal  0              0     0    0
          0             0             Geral.mtotal   0     0    0
          0             0             0              Ix   -Ixy  -Ixz
          0             0             0             -Iyx   Iy   -Iyz
          0             0             0             -Izx  -Izy   Iz  ];

    % Matriz de Coriolis
    C = [ 0            0            0                    0               Geral.mtotal*vz     -Geral.mtotal*vy
          0            0            0             -Geral.mtotal*vz             0              Geral.mtotal*vx
          0            0            0              Geral.mtotal*vy      -Geral.mtotal*vx            0
          0            0            0                    0                     0                    0
          0            0            0                 -Iy*wz                   0                  Ix*wx
          0            0            0                  Iy*wy                -Ix*wx                  0          ];

    % Matriz de forças aerodinâmicas
    D = 0.5 * rho * Corpo.AreaFoguete * [ Cd0(i+1)*Va    (Cdi(i+1)-Clalpha(i+1))*vy           (Cdi(i+1)-Clalpha(i+1))*vz          0  0  0
                                          Cd0(i+1)*vy    Clalpha(i+1)*Va+(Cdi(i+1)*vy^2/Va)   Cdi(i+1)*vz*vy/Va                   0  0  0
                                          Cd0(i+1)*vz    Cdi(i+1)*vy*vz/Va                    Clalpha(i+1)*Va+(Cdi(i+1)*vz^2/Va)  0  0  0
                                          0              0                                    0                                   0  0  0
                                          0              0                                   -Cmalpha(i+1)*cref*Va                0  0  0
                                          0              Cmalpha(i+1)*cref*Va                 0                                   0  0  0 ];

    % Matriz de força gravitacional
    G = -Geral.mtotal * g * [    -sin(phi_y)
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
iimax = 5000;   % Valor máximo da variável de controle
t_p1 = 0;       % Condição inicial de tempo do paraquedas inicial

% Condições iniciais de velocidade (referencial inercial)
eta_pontoX = -5;
eta_pontoY = eta_ponto(2,end);
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
        t(i+1) = t(i) + dt;
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
                    0   ] + cumtrapz(t_p1, eta_ponto_p1, 2);

        %% Incremento das variáveis de controle
        i = i + 1;
        ii = ii + 1;
    end
end

%% Acionamento do paraquedas principal

iii = 1;         % Variável de controle do paraquedas inicial
iiimax = 5000;   % Valor máximo da variável de controle
t_p2 = 0;        % Condição inicial de tempo do paraquedas inicial

% Condições iniciais de velocidade (referencial inercial)
eta_pontoX = -5;
eta_pontoY = eta_ponto_p1(2,end);
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
        t(i+1) = t(i) + dt;
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
        i = i + 1;
        iii = iii + 1;
    end
end
return
%% Plots

% Deslocamento 3D
plotaPos3D(t, eta, eta_p1, eta_p2)

% Módulo da velocidade
plotaModVel(t, eta_ponto, eta_ponto_p1, eta_ponto_p2, MACH)

% Ângulos de ataque e escorregamento
plotaAngulos(t, alpha_deg, beta_deg)

% Coeficientes aerodinâmico
plotaCoef(t, Cd0, Cdi, Clalpha, Cmalpha)

%% Animações

% Deslocamento 3D
animaPos3D(t, eta, eta_p1, eta_p2, eta_ponto, eta_ponto_p1, eta_ponto_p2)
