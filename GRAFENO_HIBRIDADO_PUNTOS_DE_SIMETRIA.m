clear all;
close all;
clc;

% Constantes del modelo
a = 2.46;
Es = -8.37;
Epx = 0;
Epy = 0;
Epz = 0;   
Vsssigma = -5.73;
Vspsigma = 6.05;
%Vpppi = 0;
Vpppi = -3.07;
Vppsigma = 5.62;

Esp2 = (1/3)*Es+(2/3)*Epx;
Esp2t = (1/3)*Es + (1/6)*Epx + (1/2)*Epy;
Esp2tt = (1/3)*Es + (1/6)*Epx + (1/2)*Epy;

V1arriba = (1/3)*Vsssigma-(2/3)*Vppsigma-(2*sqrt(2)/3)*Vspsigma;
V1abajo = (1/3)*Vsssigma-(2/3)*Vppsigma+(2*sqrt(2)/3)*Vspsigma;

V2 = (1/3)*Vsssigma - (1/sqrt(18))*((1/2)*Vspsigma) - (1/sqrt(6))*(sqrt(3)/2)*Vspsigma + ...
(1/sqrt(18))*(-(1/2)*Vspsigma) - (1/6)*((Vppsigma/4)+(3*Vpppi/4)) - (1/sqrt(12))*((sqrt(3)/4)*(Vppsigma-Vpppi))+ ...
(1/sqrt(6))*(-(sqrt(3)/2)*Vspsigma) - (1/sqrt(12))*((sqrt(3)/4)*(Vppsigma-Vpppi)) - (1/2)*((3*Vppsigma/4)+(Vpppi/4));

V3 = (1/3)*Vsssigma - (1/sqrt(18))*(Vspsigma/2) + (1/sqrt(6))*(-sqrt(3)/2)*Vspsigma + ...
(1/sqrt(18))*(-(1/2)*Vspsigma) - (1/6)*((Vppsigma/4)+(3*Vpppi/4)) + (1/sqrt(12))*((-sqrt(3)/4)*(Vppsigma-Vpppi)) + ...
-(1/sqrt(6))*((sqrt(3)/2)*Vspsigma) + (1/sqrt(12))*((-sqrt(3)/4)*(Vppsigma-Vpppi)) - (1/2)*((3*Vppsigma/4)+(Vpppi/4));


% Número de puntos por tramo
N = 1000;
% Puntos de alta simetría

K = (2*pi/3)*[sqrt(3), 1];
G = [0, 0];
M = (2*pi/3)*[sqrt(3), 0];

% Definir segmentos del camino de simetría
k_GK  = [linspace(G(1), K(1), N )', linspace(G(2), K(2), N )'];
k_KM  = [linspace(K(1), M(1), N )', linspace(K(2), M(2), N )'];
k_MG  = [linspace(M(1), G(1), N)', linspace(M(2), G(2), N )'];

% Concatenar camino completo
k_path = [k_GK; k_KM; k_MG];

% Prealocar bandas
bands = zeros(size(k_path,1), 8);

% Cálculo de las bandas
for i = 1:size(k_path,1)
    kx = k_path(i,1);
    ky = k_path(i,2);
   
    % Matriz h1
        h1=[Esp2 0 0 0 V1arriba 0 0 0; 
        0 Esp2t 0 0 0 0 0 0; 
        0 0 Esp2tt 0 0 0 0 0; 
        0 0 0 Epz 0 0 0 Vpppi; 
        V1abajo 0 0 0 Esp2 0 0 0; 
        0 0 0 0 0 Esp2t 0 0; 
        0 0 0 0 0 0 Esp2tt 0 ; 
        0 0 0 Vpppi 0 0 0 Epz];

        % Matriz h2
        A = zeros(8, 8);
        B = [ 0 0 0 0;
             0 0 0 0;
             0 0 V2 0;
             0 0 0 Vpppi];
        A(5:8, 1:4) = B;
        h2 = A * exp(1i * ((sqrt(3)/2)*kx + (1/2)*ky));

        % Matriz h3
        A = zeros(8, 8);
        C = [ 0 0 0 0;
             0 V3 0 0;
             0 0 0 0;
             0 0 0 Vpppi];
        A(5:8, 1:4) = C;
        h3 = A * exp(1i * ((sqrt(3)/2)*kx - (1/2)*ky));

        % Matriz h4
        h4 = h2';

        % Matriz h5
        h5 = h3';

    % Hamiltoniano total
    H = h1 + h2 + h3 + h4 + h5;
    % Calcular autovalores ordenados
    bands(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end

% Eje x para graficar
x = 1:size(k_path,1);

% Posiciones y etiquetas de los puntos de alta simetría
xticks_pos = [1, N, 2*N, 3*N];
xtick_labels = {'\Gamma','K','M','\Gamma'};

% Gráfico
figure;
plot(x, bands, 'k', 'LineWidth', 1.2); % todas líneas negras
xline(xticks_pos, '--k'); % líneas verticales en puntos de simetría
xticks(xticks_pos);
xticklabels(xtick_labels);
xlabel('\sl Trayectoria en la zona de Brillouin','FontSize', 10,'Color', 'k');
ylabel('\sl Energía (eV)','FontSize', 10,'Color', 'k');
title('');
ylim([-25 10]); % ajustar al rango típico
grid on;