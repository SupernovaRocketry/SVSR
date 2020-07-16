function animaPos3D(t, eta, eta_p1, eta_p2, eta_ponto, eta_ponto_p1, eta_ponto_p2)
%% Simplificação das coordenadas
% Posição
etaTotal = cat(2, eta, eta_p1(:,2:end), eta_p2(:,2:end));
x = etaTotal(1,:);
y = etaTotal(2,:);
z = abs(etaTotal(3,:));

% Velocidade
eta_pontoTotal = cat(2, eta_ponto, eta_ponto_p1(:,2:end), eta_ponto_p2(:,2:end));
vx = eta_pontoTotal(1,:);
vy = eta_pontoTotal(2,:);
vz = eta_pontoTotal(3,:);
V = sqrt(vx.^2 + vy.^2 + vz.^2);

%% Animação
figure
set(gcf,'Units','centimeters','Position',[8 8 30 16])
set(gcf,'color','w')
box on
dx = x(1);
dy = y(1);
dz = z(1);
for i = 2:10:size(eta,2)
    dx = [dx , x(i)];
    dy = [dy , y(i)];
    dz = [dz , z(i)];
    plot3(dx , dy , dz , 'Linewidth' , 2)
    hold on
    scatter3(x(i) , y(i) , z(i) , 'ro' , 'linewidth' , 1)
    title(sprintf('Tempo de simulação: %.0f s' , t(i)))
    zlim([0 2000])
    ylim([-100 100])
    xlim([-1000 1000])
    aa = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    bb = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    cc = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    dd = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    legend([aa bb cc dd] , sprintf('d_x = %.0f m' , x(i)) , sprintf('d_y = %.1f m' , y(i)) , sprintf('d_z = %.0f m' , z(i)) , sprintf('v = %.0f m/s' , V(i)) , 'location' , 'northwest')
    grid on
    view(-20, 30)
    drawnow
    hold off
end
