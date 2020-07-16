% C�digo modificado para c�lculo da EDO do
% modelo determin�stico do simulador com 6DOF
% pelo m�todo de Runge-Kutta de 4� ordem
%%
imax = 20;          % Valor m�ximo da vari�vel de controle
h = 0.01;           % Passo da simula��o
x(1) = 0;           % Condi��o inicial de tempo
v(:,1) = v;         % Condi��o inicial de velocidade

% Interpola��o do expuxo
Thrust = [interp1(Thrust_t , Thrust_N , 0:h:Thrust_t(end)) , zeros(1,1e6)];

% EDO do sistema
F_xy = @(t,r,tau) -(C/M) * r - (D/M) * r - (M\G) + (M\tau);

%% For�as atuantes durante lan�mento
for i = 1:10
    % Matrizes de for�as
    G = zeros(6,1);
    D = zeros(6,6);
    tauSim = [ Thrust(i+1)
                  0
                  0
                  0
                  0
                  0        ];
    %% M�todo de Runge-Kutta
    % Coeficientes do equacionamento de Runge-Kutta
    k_1 = F_xy(x(i) , v(:,i) , tauSim);
    k_2 = F_xy(x(i)+0.5*h , v(:,i)+0.5*h*k_1 , tauSim);
    k_3 = F_xy((x(i)+0.5*h) , (v(:,i)+0.5*h*k_2) , tauSim);
    k_4 = F_xy((x(i)+h) , (v(:,i)+k_3*h) , tauSim);
    % Equa��o principal de Runge-Kutta
    v(:,i+1) = v(:,i) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h;
    
    %% Atualiza��o das vari�veis da simula��o
    % Incremento do tempo
    x(i+1) = x(i) + h;

    % Vetor velocidade no triedro inercial
    eta_ponto(:,i+1) = J * v(:,i+1);

    % Vetor posi��o no triedro inercial
    eta = cumtrapz(x, eta_ponto, 2);

    % Simplifica��es de coordenadas
    vx = v(1,i+1);
    vy = v(2,i+1);
    vz = v(3,i+1);
    wx = v(4,i+1);
    wy = v(5,i+1);
    wz = v(6,i+1);

    phi_x = eta(4,i+1);
    phi_y = deg2rad(-85);
    phi_z = eta(6,i+1);

    % M�dulo da velocidade
    Va = calculaModulo(vx,vy,vz);
    % N�mero de MACH
    MACH = Va / 360;

    % �ngulos de ataque e escorregamento
    [alpha_deg(i), alpha_rad] = calculaAlpha(vz,vx);
    [beta_deg(i), beta_rad] = calculaBeta(vy,Va);

    % Coeficiente de Sustenta��o (Clalpha)
    Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans;

    % Coeficientes de Arrasto (Cd0 e Cdi)
    Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha);
    Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad, MACH, Clalpha);

    % Coeficiente de Momento Aerodin�mico (Cmalpha)
    cref = Aleta.Lr;
    Cmalpha = calculaCoefMomento(Geral, Clalpha, cref);
    
    %% Atualiza��o das matrizes
    % Matriz de transforma��o de coordenadas (Corpo > Inercial)
    J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
          sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
         -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
          0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
          0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
          0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

    % Matriz de momentos
    M = [ Geral.m  0        0         0     0    0
          0        Geral.m  0         0     0    0
          0        0        Geral.m   0     0    0
          0        0        0         Ix   -Ixy  -Ixz
          0        0        0        -Iyx   Iy   -Iyz
          0        0        0        -Izx  -Izy   Iz  ];

    % Matriz de Coriolis
    C = [ 0            0            0                   0               Geral.m*vz       -Geral.m*vy
          0            0            0              -Geral.m*vz              0             Geral.m*vx
          0            0            0               Geral.m*vy         -Geral.m*vx             0
          0            Geral.m*vz  -Geral.m*vy   Ixy*vz - Ixz*wy     Ixz*wx + Ix*wz    -Ixy*wx - Ix*wy
         -Geral.m*vz   0            Geral.m*vx   -Iy*vz - Iyz*wy     Iyz*wx - Iyx*wz     Iy*wx + Iyx*wy
          Geral.m*vy  -Geral.m*vx   0            Izy*vz + Iz*wy      -Iz*wx - Izx*wz    Izx*wy - Izy*wx ];
end
i = i + 1;

%% For�as atuantes no espa�o livre
while true
    % Defini��o da matriz tau
    tauSim = [ Thrust(i)
                  0
                  0
                  0
                  0
                  0      ];
        %% M�todo de Runge-Kutta
        if (v(3,i) < -0.1 || i == imax)  % Se vz < -0.1 (descida) ou i = imax (overflow)
            break;
        else
            % Coeficientes do equacionamento de Runge-Kutta
            k_1 = F_xy(x(i) , v(:,i) , tauSim);
            k_2 = F_xy(x(i)+0.5*h , v(:,i)+0.5*h*k_1 , tauSim);
            k_3 = F_xy((x(i)+0.5*h) , (v(:,i)+0.5*h*k_2) , tauSim);
            k_4 = F_xy((x(i)+h) , (v(:,i)+k_3*h) , tauSim);
            % Equa��o principal de Runge-Kutta
            v(:,i+1) = v(:,i) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h;
        end
        %% Atualiza��o das vari�veis da simula��o
        % Incremento do tempo
        x(i+1) = x(i) + h;

        % Vetor velocidade no triedro inercial
        eta_ponto(:,i+1) = J * v(:,i+1);

        % Vetor posi��o no triedro inercial
        eta = cumtrapz(x, eta_ponto, 2);

        % Simplifica��es de coordenadas
        vx = v(1,i+1);
        vy = v(2,i+1);
        vz = v(3,i+1);
        wx = v(4,i+1);
        wy = v(5,i+1);
        wz = v(6,i+1);

        phi_x = eta(4,i+1);
        phi_y = eta(5,i+1);
        phi_z = eta(6,i+1);

        % M�dulo da velocidade
        Va = calculaModulo(vx,vy,vz);
        % N�mero de MACH
        MACH = Va / 360;

        % �ngulos de ataque e escorregamento
        [alpha_deg(i), alpha_rad] = calculaAlpha(vz,vx);
        [beta_deg(i), beta_rad] = calculaBeta(vy,Va);

        % Coeficiente de Sustenta��o (Clalpha)
        Clalpha = Cl_coifa + Cl_corpo*alpha_rad + Cl_aleta + Cl_trans;

        % Coeficientes de Arrasto (Cd0 e Cdi)
        Cd0 = calculaArrastoMinimo(Cd_corpo, Cd_base, Cd_aleta, Cd_interf, Va, Geral, Aleta, alpha_rad, MACH, Clalpha);
        Cdi = calculaArrastoAoA(Geral, Coifa, Corpo, Aleta, alpha_rad, MACH, Clalpha);

        % Coeficiente de Momento Aerodin�mico (Cmalpha)
        cref = Aleta.Lr;
        Cmalpha = calculaCoefMomento(Geral, Clalpha, cref);

        %% Atualiza��o das matrizes
        % Matriz de transforma��o de coordenadas (Corpo > Inercial)
        J = [ cos(phi_z)*cos(phi_y) -sin(phi_z)*cos(phi_x)+cos(phi_z)*sin(phi_y)*sin(phi_x)   sin(phi_z)*sin(phi_x)+cos(phi_z)*cos(phi_x)*sin(phi_y)  0    0                      0
              sin(phi_z)*cos(phi_y)  cos(phi_z)*cos(phi_x)+sin(phi_x)*sin(phi_y)*sin(phi_z)  -cos(phi_z)*sin(phi_x)+sin(phi_y)*sin(phi_z)*cos(phi_x)  0    0                      0
             -sin(phi_y)             cos(phi_y)*sin(phi_x)                                    cos(phi_y)*cos(phi_x)                                   0    0                      0
              0                      0                                                        0                                                       1    sin(phi_x)*tan(phi_y)  cos(phi_x)*tan(phi_y)
              0                      0                                                        0                                                       0    cos(phi_x)            -sin(phi_x)
              0                      0                                                        0                                                       0    sin(phi_x)/cos(phi_y)  cos(phi_x)/cos(phi_y) ];

        % Matriz de momentos
        M = [ Geral.m  0        0         0     0    0
              0        Geral.m  0         0     0    0
              0        0        Geral.m   0     0    0
              0        0        0         Ix   -Ixy  -Ixz
              0        0        0        -Iyx   Iy   -Iyz
              0        0        0        -Izx  -Izy   Iz  ];

        % Matriz de Coriolis
        C = [ 0            0            0                   0               Geral.m*vz       -Geral.m*vy
              0            0            0              -Geral.m*vz              0             Geral.m*vx
              0            0            0               Geral.m*vy         -Geral.m*vx             0
              0            Geral.m*vz  -Geral.m*vy   Ixy*vz - Ixz*wy     Ixz*wx + Ix*wz    -Ixy*wx - Ix*wy
             -Geral.m*vz   0            Geral.m*vx   -Iy*vz - Iyz*wy     Iyz*wx - Iyx*wz     Iy*wx + Iyx*wy
              Geral.m*vy  -Geral.m*vx   0            Izy*vz + Iz*wy      -Iz*wx - Izx*wz    Izx*wy - Izy*wx ];

        % Matriz de for�as aerodin�micas
        D = 0.5 * rho * Corpo.AreaFoguete * [ Cd0*Va   (Cdi-Clalpha)*vy          (Cdi-Clalpha)*vz          0  0  0  0  0  0
                                              Cd0*vy   Clalpha*Va+(Cdi*vy^2/Va)  Cdi*vz*vy/Va              0  0  0  0  0  0
                                              Cd0*vz   Cdi*vy*vz/Va              Clalpha*Va+(Cdi*vz^2/Va)  0  0  0  0  0  0
                                              0        0                         0                         0  0  0  0  0  0
                                              0        0                        -Cmalpha*cref*Va           0  0  0  0  0  0
                                              0        Cmalpha*cref*Va           0                         0  0  0  0  0  0 ];

        % Matriz de for�a gravitacional
        G = -Geral.m * g * [    -sin(phi_x)
                            cos(phi_y)*sin(phi_x)
                            cos(phi_y)*cos(phi_x)
                                     0
                                     0
                                     0            ];

    %% Incremento da vari�vel de controle
    i = i + 1;
end
