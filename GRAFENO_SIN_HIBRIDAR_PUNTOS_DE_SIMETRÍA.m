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
Vpppi = -3.07;
Vppsigma = 5.62;

% Número de puntos por tramo
N = 100;
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
        h1=[Es 0 0 0 Vsssigma Vspsigma 0 0; 
        0 Epx 0 0 -Vspsigma Vppsigma 0 0; 
        0 0 Epy 0 0 0 Vpppi 0; 
        0 0 0 Epz 0 0 0 Vpppi; 
        Vsssigma -Vspsigma 0 0 Es 0 0 0; 
        Vspsigma Vppsigma 0 0 0 Epx 0 0; 
        0 0 Vpppi 0 0 0 Epy 0 ; 
        0 0 0 Vpppi 0 0 0 Epz];

        % Matriz h2
        A = zeros(8, 8);
        B = [ Vsssigma (1/2)*Vspsigma (sqrt(3)/2)*Vspsigma 0;
             -(1/2)*Vspsigma (Vppsigma/4)+(3*Vpppi/4) (sqrt(3)/4)*(Vppsigma-Vpppi) 0;
             -(sqrt(3)/2)*Vspsigma (sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0;
             0 0 0 Vpppi];
        A(5:8, 1:4) = B;
        h2 = A * exp(1i * ((sqrt(3)/2)*kx + (1/2)*ky));

        % Matriz h3
        A = zeros(8, 8);
        C = [
        Vsssigma (Vspsigma/2) -(sqrt(3)/2)*Vspsigma 0 ; 
        -(Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) -(sqrt(3)/4)*(Vppsigma-Vpppi) 0 ;
        (sqrt(3)/2)*Vspsigma -(sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0 ; 
        0 0 0 Vpppi]; 
        A(5:8, 1:4) = C;
        h3 = A * exp(1i * ((sqrt(3)/2)*kx - (1/2)*ky));

        % Matriz h4
        A = zeros(8, 8);
        D = [
        Vsssigma -(Vspsigma/2) -(sqrt(3)/2)*Vspsigma 0; 
        (Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) (sqrt(3)/4)*(Vppsigma-Vpppi) 0; 
        (sqrt(3)/2)*Vspsigma (sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0;
        0 0 0 Vpppi];
        A(1:4, 5:8) = D;
        h4 = A * exp(1i * ((-sqrt(3)/2)*kx - (1/2)*ky));

        % Matriz h5
        A = zeros(8, 8);
        E = [Vsssigma -(Vspsigma/2) (sqrt(3)/2)*Vspsigma 0; 
        (Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) -(sqrt(3)/4)*(Vppsigma-Vpppi) 0; 
        -(sqrt(3)/2)*Vspsigma -(sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0; 
        0 0 0 Vpppi];
        A(1:4, 5:8) = E;
        h5 = A * exp(1i * ((-sqrt(3)/2)*kx + (1/2)*ky));

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
ylim([-28 15]); % ajustar al rango típico
grid on;