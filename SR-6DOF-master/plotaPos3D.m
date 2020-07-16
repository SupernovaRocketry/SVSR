function plotaPos3D(t, eta, eta_p1, eta_p2)
%% Simplificação das coordenadas
etaTotal = cat(2, eta, eta_p1(:,2:end), eta_p2(:,2:end));

x = etaTotal(1,:);
y = etaTotal(2,:);
z = abs(etaTotal(3,:));

%% Plot Total
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot3(x , y , z , 'Linewidth' , 2)
zlim([0 2000])
ylim([-100 100])
xlim([-1000 1000])
grid on
box on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view(-20, 30)

%% Plot por Fases de Vôo
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot3(eta(1,:) , eta(2,:) , abs(eta(3,:)) , 'color' , [0 0.447 0.741] , 'Linewidth' , 2)
hold on
plot3(eta_p1(1,2:end) , eta_p1(2,2:end) , abs(eta_p1(3,2:end)) , 'r' , 'Linewidth' , 2)
plot3(eta_p2(1,2:end) , eta_p2(2,2:end) , abs(eta_p2(3,2:end)) , 'color' , [0 0.7 0] , 'Linewidth' , 2)
hold off
zlim([0 2000])
ylim([-100 100])
xlim([-1000 1000])
grid on
box on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
view(-20, 30)
legend('Liftoff + Coasting','First Deployment Event','Second Deployment Event')
