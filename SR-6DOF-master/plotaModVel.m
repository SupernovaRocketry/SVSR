function plotaModVel(t, eta_ponto, eta_ponto_p1, eta_ponto_p2, MACH)
%% Simplificação das coordenadas
eta_pontoTotal = cat(2, eta_ponto, eta_ponto_p1(:,2:end), eta_ponto_p2(:,2:end));

vx = eta_pontoTotal(1,:);
vy = eta_pontoTotal(2,:);
vz = eta_pontoTotal(3,:);

V = sqrt(vx.^2 + vy.^2 + vz.^2);

%% Velocidade
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t , V , 'Linewidth' , 1.5)
grid on
ylabel('Velocidade (m/s)')
xlabel('Tempo (s)')
%% MACH
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(MACH)) , MACH , 'Linewidth' , 1.5)
grid on
ylabel('MACH')
xlabel('Tempo (s)')
