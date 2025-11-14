clear all;
close all;
clc;
 
% Parámetros de la red
a = 4.62;
b = 3.30;
gamma = 76.31;
alfa = 98.15/2;
s = sind(gamma);
c = cosd(gamma);
smedios = sind(alfa);
cmedios = cosd(alfa);
Omega = sqrt(1-s^2);
Delta = (-s);
Es = -8.8;
Epx = 0;
Epy = 0;
Epz = 0;   
Vsssigma = -1.59;
Vspsigma = 2.39;
Vppsigma = 4.03;
Vpppi = -1.14;

% Número de puntos por tramo
N = 1000;

% Puntos de alta simetría
G = [0, 0];
X = [pi/a, 0];
S = [pi/a, pi/b];
Y = [0, pi/b];

% Definir segmentos del camino de simetría
k_GX  = [linspace(G(1), X(1), N )', linspace(G(2), X(2), N )'];
k_XS  = [linspace(X(1), S(1), N )', linspace(X(2), S(2), N )'];
k_SY  = [linspace(S(1), Y(1), N)', linspace(S(2), Y(2), N )'];
k_YG  = [linspace(Y(1), G(1), N)', linspace(Y(2), G(2), N )'];

% Camino completo
k_path = [k_GX; k_XS; k_SY; k_YG];

% Cálculo de las bandas
bands = zeros(size(k_path,1), 16);
for i = 1:size(k_path,1)
    kx = k_path(i,1);
    ky = k_path(i,2);
   
     % Matriz h0
        h0=[
        Es 0 0 0 Vsssigma Omega*Vspsigma 0 -Delta*Vspsigma 0 0 0 0 Vsssigma -cmedios*Vspsigma -smedios*Vspsigma 0; 
        0 Epx 0 0 -Omega*Vspsigma Omega^2*Vppsigma+(1-Omega^2)*Vpppi 0 -Omega*Delta*(Vppsigma-Vpppi) 0 0 0 0 cmedios*Vspsigma cmedios^2*Vppsigma+(1-cmedios^2)*Vpppi smedios*cmedios*(Vppsigma-Vpppi) 0;
        0 0 Epy 0 0 0 Vpppi 0 0 0 0 0 smedios*Vspsigma smedios*cmedios*(Vppsigma-Vpppi) smedios^2*Vppsigma+(1-smedios^2)*Vpppi 0; 
        0 0 0 Epz Delta*Vspsigma -Omega*Delta*(Vppsigma-Vpppi) 0 Delta^2*Vppsigma+(1-Delta^2)*Vpppi 0 0 0 0 0 0 0 Vpppi; 
        Vsssigma -Omega*Vspsigma 0 Delta*Vspsigma Es 0 0 0 Vsssigma cmedios*Vspsigma -smedios*Vspsigma 0 0 0 0 0; 
        Omega*Vspsigma Omega^2*Vppsigma+(1-Omega^2)*Vpppi 0 -Omega*Delta*(Vppsigma-Vpppi) 0 Epx 0 0 -cmedios*Vspsigma cmedios^2*Vppsigma+(1-cmedios^2)*Vpppi -smedios*cmedios*(Vppsigma-Vpppi) 0 0 0 0 0; 
        0 0 Vpppi 0 0 0 Epy 0 smedios*Vspsigma -smedios*cmedios*(Vppsigma-Vpppi) smedios^2*Vppsigma+(1-smedios^2)*Vpppi 0 0 0 0 0;
        -Delta*Vspsigma -Omega*Delta*(Vppsigma-Vpppi) 0 Delta^2*Vppsigma+(1-Delta^2)*Vpppi 0 0 0 Epz 0 0 0 Vpppi 0 0 0 0; 
        0 0 0 0 Vsssigma -cmedios*Vspsigma smedios*Vspsigma 0 Es 0 0 0 0 0 0 0; 
        0 0 0 0 cmedios*Vspsigma cmedios^2*Vppsigma+(1-cmedios^2)*Vpppi -smedios*cmedios*(Vppsigma-Vpppi) 0 0 Epx 0 0 0 0 0 0; 
        0 0 0 0 -smedios*Vspsigma -smedios*cmedios*(Vppsigma-Vpppi) smedios^2*Vppsigma+(1-smedios^2)*Vpppi 0 0 0 Epy 0 0 0 0 0;  
        0 0 0 0 0 0 0 Vpppi 0 0 0 Epz 0 0 0 0; 
        Vsssigma cmedios*Vspsigma smedios*Vspsigma 0 0 0 0 0 0 0 0 0 Es 0 0 0; 
        -cmedios*Vspsigma cmedios^2*Vppsigma+(1-cmedios^2)*Vpppi smedios*cmedios*(Vppsigma-Vpppi) 0 0 0 0 0 0 0 0 0 0 Epx 0 0;
        -smedios*Vspsigma smedios*cmedios*(Vppsigma-Vpppi) smedios^2*Vppsigma+(1-smedios^2)*Vpppi 0 0 0 0 0 0 0 0 0 0 0 Epy 0; 
        0 0 0 Vpppi 0 0 0 0 0 0 0 0 0 0 0 Epz; 
        ];
        % Matriz h1
        A=[
        0 0 0 0 0 0 0 0 0 0 0 0 Vsssigma -cmedios*Vspsigma smedios*Vspsigma 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 cmedios*Vspsigma cmedios^2*Vppsigma+(1-cmedios^2)*Vpppi -smedios*cmedios*(Vppsigma-Vpppi) 0;
        0 0 0 0 0 0 0 0 0 0 0 0 -smedios*Vspsigma -smedios*cmedios*(Vppsigma-Vpppi) smedios^2*Vppsigma+(1-smedios^2)*Vpppi 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 Vpppi; 
        0 0 0 0 0 0 0 0 Vsssigma cmedios*Vspsigma smedios*Vspsigma 0 0 0 0 0; 
        0 0 0 0 0 0 0 0 -cmedios*Vspsigma cmedios^2*Vppsigma+(1-cmedios^2)*Vpppi smedios*cmedios*(Vppsigma-Vpppi) 0 0 0 0 0;
        0 0 0 0 0 0 0 0 -smedios*Vspsigma smedios*cmedios*(Vppsigma-Vpppi) smedios^2*Vppsigma+(1-smedios^2)*Vpppi 0 0 0 0 0; 
        0 0 0 0 0 0 0 0 0 0 0 Vpppi 0 0 0 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        ];
        h1 = A * exp(1i*b*ky);
        % Matriz h2
        B = [
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 Vsssigma Omega*Vspsigma 0 Delta*Vspsigma; 
        0 0 0 0 0 0 0 0 0 0 0 0 -Omega*Vspsigma Omega^2*Vppsigma+(1-Omega^2)*Vpppi 0 Delta*Omega*(Vppsigma-Vpppi);
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 Vpppi 0; 
        0 0 0 0 0 0 0 0 0 0 0 0 -Delta*Vspsigma Delta*Omega*(Vppsigma-Vpppi) 0 Delta^2*Vppsigma+(1-Delta^2)*Vpppi; 
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
        ];
        h2 = B*exp(1i*a*kx);
      
        % Hamiltoniano total
        H = h0 + h1 + h2 + h1' + h2';
    % Calcular autovalores ordenados
    bands(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end
% Eje x para graficar
x = 1:size(k_path,1);
% Posiciones y etiquetas de los puntos de alta simetría
xticks_pos = [1, N, 2*N, 3*N, 4*N];
xtick_labels = {'\Gamma','X','S','Y', '\Gamma'};
% Gráfico
figure;
plot(x, bands, 'k', 'LineWidth', 1.2); % todas líneas negras
xline(xticks_pos, '--k'); % líneas verticales en puntos de simetría
xticks(xticks_pos);
xticklabels(xtick_labels);
xlabel('\sl Trayectoria en la zona de Brillouin','FontSize', 10,'Color', 'k');
ylabel('\sl Energía (eV)','FontSize', 10,'Color', 'k');
title('');
ylim([-20 10]); % ajustar al rango típico
grid on;