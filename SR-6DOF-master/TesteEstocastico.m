close all
clear all
clc
tic
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

%% *******************************************************
%% Modelo Estocástico
%% *******************************************************

% Número de simulações
n = 2;

for k = 1:n

Geral.mtot = Geral.mtotal + (sqrt(0.1)*randn);

%% Momentos de Inércia

Ix = 687989740.24e-9 + (sqrt(0.01)*randn);
Iy = 812188074.87e-9 + (sqrt(0.01)*randn);
Iz = 126562111.96e-9 + (sqrt(0.01)*randn);
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
Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans + (sqrt(5)*randn);

% Coeficiente de Arrasto (Cd0)
Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha) + (sqrt(0.01)*randn);
Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad) + (sqrt(0.1)*randn);

% Coeficiente de Momento Aerodinâmico (Cmalpha)
cref = Aleta.Lr;
Cmalpha = calculaCoefMomento(Geral, Clalpha, cref) + (sqrt(5)*randn);

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
eta_ponto(:,1,k) = J * v_i;

%% *******************************************************
%% Resolução do Método
%% *******************************************************

imax = 5000;         % Valor máximo da variável de controle
dt = 0.01;           % Passo da simulação
t(1) = 0;            % Condição inicial de tempo
v(:,1,k) = v_i;      % Condição inicial de velocidade

% Interpolação do expuxo
Thrust = [interp1(Prop.Thrust_t , Prop.Thrust_N , 0:dt:Prop.Thrust_t(end)) , zeros(1,1e6)];
Thrust = Thrust + (sqrt(200)*randn(1,length(Thrust)));

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
    k_1 = F_xy(t(i) , v(:,i,k) , tauSim);
    k_2 = F_xy(t(i)+0.5*dt , v(:,i,k)+0.5*dt*k_1 , tauSim);
    k_3 = F_xy((t(i)+0.5*dt) , (v(:,i,k)+0.5*dt*k_2) , tauSim);
    k_4 = F_xy((t(i)+dt) , (v(:,i,k)+k_3*dt) , tauSim);
    % Equação principal de Runge-Kutta
    v(:,i+1,k) = v(:,i,k) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*dt;
    
    %% Atualização das variáveis da simulação
    % Incremento do tempo
    t(i+1) = t(i) + dt;

    % Condição da base de lançamento (velocidades perpendiculares = 0 ; velocidades angulares = 0)
    v(2,i+1,k) = 0;
    v(3,i+1,k) = 0;
    v(4,i+1,k) = 0;
    v(5,i+1,k) = 0;
    v(6,i+1,k) = 0;

    % Vetor velocidade no triedro inercial
    eta_ponto(:,i+1,k) = J * v(:,i+1,k);

    % Vetor posição no triedro inercial
    eta_int = cumtrapz(t , eta_ponto(:,:,k), 2);
    eta(:,i+1,k) = eta_int(:,i+1) + eta_i;

    % Simplificações de coordenadas
    vx = v(1,i+1,k);
    vy = v(2,i+1,k);
    vz = v(3,i+1,k);
    wx = v(4,i+1,k);
    wy = v(5,i+1,k);
    wz = v(6,i+1,k);

    phi_x = eta(4,i+1,k);
    phi_y = deg2rad(Langle);
    phi_z = eta(6,i+1,k);

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
    if (eta_ponto(3,i,k) > 0.1 || i == imax)  % Se vz > 0.1 (descida) ou i = imax (overflow)
        break;
    else
        % Coeficientes do equacionamento de Runge-Kutta
        k_1 = F_xy(t(i) , v(:,i,k) , tauSim);
        k_2 = F_xy(t(i)+0.5*dt , v(:,i,k)+0.5*dt*k_1 , tauSim);
        k_3 = F_xy((t(i)+0.5*dt) , (v(:,i,k)+0.5*dt*k_2) , tauSim);
        k_4 = F_xy((t(i)+dt) , (v(:,i,k)+k_3*dt) , tauSim);
        % Equação principal de Runge-Kutta
        v(:,i+1,k) = v(:,i,k) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*dt;
    end
    %% Atualização das variáveis da simulação
    % Incremento do tempo
    t(i+1) = t(i) + dt;

    % Vetor velocidade no triedro inercial
    eta_ponto(:,i+1,k) = J * v(:,i+1,k);

    % Vetor posição no triedro inercial
    eta_int = cumtrapz(t, eta_ponto(:,:,k), 2);
    eta(:,i+1,k) = eta_int(:,i+1) + eta_i;

    % Simplificações de coordenadas
    vx = v(1,i+1,k);
    vy = v(2,i+1,k);
    vz = v(3,i+1,k);
    wx = v(4,i+1,k);
    wy = v(5,i+1,k);
    wz = v(6,i+1,k);

    phi_x = eta(4,i+1,k);
    phi_y = eta(5,i+1,k);
    phi_z = eta(6,i+1,k);

    % Módulo da velocidade
    Va = calculaModulo(vx,vy,vz);
    % Número de MACH
    MACH = Va / 360;

    % Ângulos de ataque e escorregamento
    [alpha_deg, alpha_rad] = calculaAlpha(vz,vx);
    [beta_deg, beta_rad] = calculaBeta(vy,Va);

    % Coeficiente de Sustentação (Clalpha)
    Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans + (sqrt(5)*randn);

    % Coeficiente de Arrasto (Cd0)
    Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha) + (sqrt(0.01)*randn);
    Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad) + (sqrt(0.1)*randn);

    % Coeficiente de Momento Aerodinâmico (Cmalpha)
    cref = Aleta.Lr;
    Cmalpha = calculaCoefMomento(Geral, Clalpha, cref) + (sqrt(5)*randn);

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
end
toc
return
%% Plots

% Deslocamento 3D
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
hold on
for i = 1:n
    plot3(eta(1,all(eta(3,:,i),1),i) , eta(2,all(eta(3,:,i),1),i) , abs(eta(3,all(eta(3,:,i),1),i)) , 'color' , [0 0.447 0.741] , 'Linewidth' , 1)
    grid on
end
box on
zlabel('z (m)')
ylabel('y (m)')
xlabel('x (m)')
view(-20, 30)

% Histograma dos apogeus
hmax = max(abs(eta(3,:,:)),[],2);

pdfh = normpdf(500:3000,mean(hmax),std(hmax));          % PDF de uma Gaussiana

c = mle(permute(hmax,[3 2 1]),'distribution','logn');   % PDF de uma
pdflogh = lognpdf(500:3000,c(1),c(2));                  % LogNormal

figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
histogram(hmax , 300 , 'normalization' , 'pdf')
hold on
p1 = plot(500:3000 , pdfh , 'Linewidth' , 2);
p2 = plot(500:3000 , pdflogh , 'Linewidth' ,2);
hold off
ylabel('PDF_{h_{max}}')
xlabel('h_{max}')
grid on
legend([p1 p2] , ['Gaussiana, \mu = ' num2str(mean(hmax)) ', \sigma = ' num2str(std(hmax))] , ['LogNormal, \mu = ' num2str(c(1)) ', \sigma = ' num2str(c(2))])

% Histograma das velocidades
V = sqrt(eta_ponto(1,:,:).^2 + eta_ponto(2,:,:).^2 + eta_ponto(3,:,:).^2);
vmax = max(V,[],2);

pdfv = normpdf(100:300,mean(vmax),std(vmax));           % PDF de uma Gaussiana

c = mle(permute(vmax,[3 2 1]),'distribution','logn');   % PDF de uma
pdflogv = lognpdf(100:300,c(1),c(2));                   % LogNormal

figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
histogram(vmax , 200 , 'normalization' , 'pdf')
hold on
plot(100:300 , pdfv , 'Linewidth' , 2)
p1 = plot(100:300 , pdfv , 'Linewidth' , 2);
p2 = plot(100:300 , pdflogv , 'Linewidth' ,2);
hold off
ylabel('PDF_{v_{max}}')
xlabel('v_{max}')
grid on
legend([p1 p2] , ['Gaussiana, \mu = ' num2str(mean(vmax)) ', \sigma = ' num2str(std(vmax))] , ['LogNormal, \mu = ' num2str(c(1)) ', \sigma = ' num2str(c(2))])
