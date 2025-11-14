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

kx=1;
ky=1;

% Puntos en el espacio k
kx_vals = linspace(-pi, pi, 1000);
ky_vals = linspace(-pi, pi, 1000);

% Inicialización de la matriz para las energías
bandas = zeros(16, length(kx_vals), length(ky_vals));

% Bucle sobre el espacio k
for ix = 1:length(kx_vals)
    for iy = 1:length(ky_vals)

        kx = kx_vals(ix);
        ky = ky_vals(iy);

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

        % Calcular autovalores
        bandas(:, ix, iy) = real(eig(H));

    end
end

%Visualizar las bandas
figure;
hold on;
for n = 1:16
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
