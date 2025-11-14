clear all;
close all;
clc;

% Parámetros
a = 2.46;
Es = -8.37;
Epx = 0;
Epy = 0;
Epz = 0;   
Vsssigma = -5.73;
Vspsigma = 6.05;
Vpppi = -3.07;
Vppsigma = 5.62;

kx=1;
ky=1;

% Puntos en el espacio k
kx_vals = linspace(-pi, pi, 1000);
ky_vals = linspace(-pi, pi, 1000);

% Inicialización de la matriz para las energías
bandas = zeros(8, length(kx_vals), length(ky_vals));

% Bucle sobre el espacio k
for ix = 1:length(kx_vals)
    for iy = 1:length(ky_vals)

        kx = kx_vals(ix);
        ky = ky_vals(iy);

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
        h2 = A * exp(a*1i * ((sqrt(3)/2)*kx + (1/2)*ky));

        % Matriz h3
        A = zeros(8, 8);
        C = [
        Vsssigma (Vspsigma/2) -(sqrt(3)/2)*Vspsigma 0 ; 
        -(Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) -(sqrt(3)/4)*(Vppsigma-Vpppi) 0 ;
        (sqrt(3)/2)*Vspsigma -(sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0 ; 
        0 0 0 Vpppi]; 
        A(5:8, 1:4) = C;
        h3 = A * exp(a*1i * ((sqrt(3)/2)*kx - (1/2)*ky));

        % Matriz h4
        A = zeros(8, 8);
        D = [
        Vsssigma -(Vspsigma/2) -(sqrt(3)/2)*Vspsigma 0; 
        (Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) (sqrt(3)/4)*(Vppsigma-Vpppi) 0; 
        (sqrt(3)/2)*Vspsigma (sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0;
        0 0 0 Vpppi];
        A(1:4, 5:8) = D;
        h4 = A * exp(a*1i * ((-sqrt(3)/2)*kx - (1/2)*ky));

        % Matriz h5
        A = zeros(8, 8);
        E = [Vsssigma -(Vspsigma/2) (sqrt(3)/2)*Vspsigma 0; 
        (Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) -(sqrt(3)/4)*(Vppsigma-Vpppi) 0; 
        -(sqrt(3)/2)*Vspsigma -(sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0; 
        0 0 0 Vpppi];
        A(1:4, 5:8) = E;
        h5 = A * exp(a*1i * ((-sqrt(3)/2)*kx + (1/2)*ky));

        % Hamiltoniano total
        H = h1 + h2 + h3 + h4 + h5;
        % Calcular autovalores
        bandas(:, ix, iy) = real(eig(H));

    end
end

%Visualizar las bandas
figure;
hold on;
for n = 1:8
    surf(kx_vals, ky_vals, squeeze(bandas(n, :, :)), 'EdgeColor', 'none');
end
xlabel('\sl k_x','FontSize', 20,'Color', 'k');
ylabel('\sl k_y','FontSize', 20,'Color', 'k');
zlabel('\sl E','FontSize', 20,'Color', 'k');
title('');
colormap jet;
colorbar;
view(3);
grid on;
hold off;
