clear all;
close all;
clc;

% Parámetros
a = 2.46;

Es = -3.204545;
Epx = 0;
Epy = 0;
Epz = 0; 
Vsssigma = -5.71;
Vspsigma = 5.42;
Vpppi = -3.07;
Vppsigma = 6.20;

Esp2 = (1/3)*Es+(2/3)*Epx;
Esp2t = (1/3)*Es + (1/6)*Epx + (1/2)*Epy;
Esp2tt = (1/3)*Es + (1/6)*Epx + (1/2)*Epy;

V1arriba = (1/3)*Vsssigma - (2/3)*Vppsigma - (2*sqrt(2)/3)*Vspsigma;
V1abajo = (1/3)*Vsssigma - (2/3)*Vppsigma - (2*sqrt(2)/3)*Vspsigma;

V2 = (1/3)*Vsssigma - (1/sqrt(18))*((1/2)*Vspsigma) - (1/sqrt(6))*(sqrt(3)/2)*Vspsigma + ...
(1/sqrt(18))*(-(1/2)*Vspsigma) - (1/6)*((Vppsigma/4)+(3*Vpppi/4)) - (1/sqrt(12))*((sqrt(3)/4)*(Vppsigma-Vpppi))+ ...
(1/sqrt(6))*(-(sqrt(3)/2)*Vspsigma) - (1/sqrt(12))*((sqrt(3)/4)*(Vppsigma-Vpppi)) - (1/2)*((3*Vppsigma/4)+(Vpppi/4));

V3 = (1/3)*Vsssigma - (1/sqrt(18))*(Vspsigma/2) + (1/sqrt(6))*(-sqrt(3)/2)*Vspsigma + ...
(1/sqrt(18))*(-(1/2)*Vspsigma) - (1/6)*((Vppsigma/4)+(3*Vpppi/4)) - (1/sqrt(12))*((-sqrt(3)/4)*(Vppsigma-Vpppi)) + ...
-(1/sqrt(6))*((sqrt(3)/2)*Vspsigma) + (1/sqrt(12))*((-sqrt(3)/4)*(Vppsigma-Vpppi)) - (1/2)*((3*Vppsigma/4)+(Vpppi/4));

kx=pi;
ky=pi;

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
        h2 = A * exp(a*1i * ((sqrt(3)/2)*kx + (1/2)*ky));

        % Matriz h3
        A = zeros(8, 8);
        C = [ 0 0 0 0;
             0 V3 0 0;
             0 0 0 0;
             0 0 0 Vpppi];
        A(5:8, 1:4) = C;
        h3 = A * exp(a*1i * ((sqrt(3)/2)*kx - (1/2)*ky));

        % Matriz h4
        h4 = h2';

        % Matriz h5
        h5 = h3';

        % Hamiltoniano total
        H = h1 + h2 + h3 + h4 + h5;

        autovalores = real(eig(H));
        disp(autovalores);

        %disp(h2-h4');
        %disp(h3-h5');
        %disp(h1);
        %disp(h1');
        %disp(h1-h1');


       A = H-H';
        disp(A);