function plotaCoef(t, Cd0, Cdi, Clalpha, Cmalpha)
%% Coeficiente de arrasto m�nimo
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(Cd0)) , Cd0 , 'Linewidth' , 1.5)
grid on
ylabel('Cd_0')
xlabel('Tempo (s)')
%% Coeficiente de arrasto de interfer�ncia
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(Cdi)) , Cdi , 'Linewidth' , 1.5)
grid on
ylabel('Cd_i')
xlabel('Tempo (s)')
%% Coeficiente de sustenta��o
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(Clalpha)) , Clalpha , 'Linewidth' , 1.5)
grid on
ylabel('C_{L\alpha}')
xlabel('Tempo (s)')
%% Coeficiente de momento aerodin�mico
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(t(1:length(Cmalpha)) , Cmalpha , 'Linewidth' , 1.5)
grid on
ylabel('C_{M\alpha}')
xlabel('Tempo (s)')
%% Curva Polar dos Coeficientes
figure
set(gcf,'color','w')
set(gcf,'Units','centimeters','Position',[8 8 25 14])
plot(Cd0+Cdi , Clalpha , 'Linewidth' , 1.5)
grid on
ylabel('C_D')
xlabel('C_L')
