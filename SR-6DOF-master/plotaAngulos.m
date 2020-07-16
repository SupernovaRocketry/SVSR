function plotaAngulos(t, alpha, beta)
%% Ângulo de ataque
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(alpha)) , alpha , 'Linewidth' , 1.5)
grid on
ylabel('AoA (graus)')
xlabel('Tempo (s)')
%% Ângulo de escorregamento
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(beta)) , beta , 'Linewidth' , 1.5)
grid on
ylabel('Escorregamento (graus)')
xlabel('Tempo (s)')
