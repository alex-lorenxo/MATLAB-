clear all;
close all;
clc;

% Parámetros
a = 2.46;
Es = -3.204545;
Epx = 2.690589;
Epy = 3.670414;
Epz = 0;   
Vsssigma = -5.71;
Vspsigma = 5.42;
Vpppi = -3.07;
Vppsigma = 6.20;
kx=pi;
ky=pi;

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
        h2 = A * exp(a*1i * ((1/2)*kx + (sqrt(3)/2)*ky));

        % Matriz h3
        A = zeros(8, 8);
        C = [
        Vsssigma (Vspsigma/2) -(sqrt(3)/2)*Vspsigma 0 ; 
        -(Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) -(sqrt(3)/4)*(Vppsigma-Vpppi) 0 ;
        (sqrt(3)/2)*Vspsigma -(sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0 ; 
        0 0 0 Vpppi]; 
        A(5:8, 1:4) = C;
        h3 = A * exp(a*1i * ((1/2)*kx - (sqrt(3)/2)*ky));

        % Matriz h4
        A = zeros(8, 8);
        D = [
        Vsssigma -(Vspsigma/2) -(sqrt(3)/2)*Vspsigma 0; 
        (Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) (sqrt(3)/4)*(Vppsigma-Vpppi) 0; 
        (sqrt(3)/2)*Vspsigma (sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0;
        0 0 0 Vpppi];
        A(1:4, 5:8) = D;
        h4 = A * exp(-a*1i * ((1/2)*kx + (sqrt(3)/2)*ky));

        % Matriz h5
        A = zeros(8, 8);
        E = [Vsssigma -(Vspsigma/2) (sqrt(3)/2)*Vspsigma 0; 
        (Vspsigma/2) (Vppsigma/4)+(3*Vpppi/4) -(sqrt(3)/4)*(Vppsigma-Vpppi) 0; 
        -(sqrt(3)/2)*Vspsigma -(sqrt(3)/4)*(Vppsigma-Vpppi) (3*Vppsigma/4)+(Vpppi/4) 0; 
        0 0 0 Vpppi];
        A(1:4, 5:8) = E;
        h5 = A * exp(-a*1i * ((1/2)*kx - (sqrt(3)/2)*ky));

        % Hamiltoniano total
        H = h1 + h2 + h3 + h4 + h5;

        autovalores = real(eig(H));
        A = H-H';
        disp(A);