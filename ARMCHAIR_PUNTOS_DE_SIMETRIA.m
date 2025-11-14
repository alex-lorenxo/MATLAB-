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
N = 1000;

% Extremos de cada segmento
primeroi = [0, -pi]; 
primerof = [0, pi];
segundoi = [(2*pi/(3*sqrt(3))), -pi]; 
segundof = [(2*pi/(3*sqrt(3))), pi];
terceroi = [(4*pi/(3*sqrt(3))), -pi]; 
tercerof = [(4*pi/(3*sqrt(3))), pi];
cuartoi = [(6*pi/(3*sqrt(3))), -pi]; 
cuartof = [(6*pi/(3*sqrt(3))), pi];
quintoi = [(8*pi/(3*sqrt(3))), -pi]; 
quintof = [(8*pi/(3*sqrt(3))), pi];
sextoi = [(10*pi/(3*sqrt(3))), -pi]; 
sextof = [(10*pi/(3*sqrt(3))), pi];

% Definir segmentos del camino de simetría

k_1 = [linspace(primeroi(1), primerof(1), N )', linspace(primeroi(2), primerof(2), N )'];
k_2 = [linspace(segundoi(1), segundof(1), N )', linspace(segundoi(2), segundof(2), N )'];
k_3 = [linspace(terceroi(1), tercerof(1), N )', linspace(terceroi(2), tercerof(2), N )'];
k_4 = [linspace(cuartoi(1), cuartof(1), N )', linspace(cuartoi(2), cuartof(2), N )'];
k_5 = [linspace(quintoi(1), quintof(1), N )', linspace(quintoi(2), quintof(2), N )'];
k_6 = [linspace(sextoi(1), sextof(1), N )', linspace(sextoi(2), sextof(2), N )'];

% Prealocar bandas
bands1 = zeros(size(k_1, 1), 8);
bands2 = zeros(size(k_2, 1), 8);
bands3 = zeros(size(k_3, 1), 8);
bands4 = zeros(size(k_4, 1), 8);
bands5 = zeros(size(k_5, 1), 8);
bands6 = zeros(size(k_6, 1), 8);

% Cálculo de las bandas
for i = 1:size(k_1,1)
    kx = k_1(i,1);
    ky = k_1(i,2);
   
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
    bands1(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end

% Cálculo de las bandas
for i = 1:size(k_2,1)
    kx = k_2(i,1);
    ky = k_2(i,2);
   
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
    bands2(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end

% Cálculo de las bandas
for i = 1:size(k_3,1)
    kx = k_3(i,1);
    ky = k_3(i,2);
   
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
    bands3(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end

% Cálculo de las bandas
for i = 1:size(k_4,1)
    kx = k_4(i,1);
    ky = k_4(i,2);
   
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
    bands4(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end

% Cálculo de las bandas
for i = 1:size(k_5,1)
    kx = k_5(i,1);
    ky = k_5(i,2);
   
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
    bands1(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end

% Cálculo de las bandas
for i = 1:size(k_6,1)
    kx = k_6(i,1);
    ky = k_6(i,2);
   
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
    bands6(i,:) = sort(real(eig(H)))'; % usar sort para continuidad visual
end


% Eje x para graficar
x = 1:size(k_1,1);

% Gráfico
figure;
plot(x, bands1, 'k', x, bands2, 'k', x, bands3, 'k', x, bands4, 'k', x, bands5, 'k', x, bands6, 'k', 'LineWidth', 1.2); hold on;


% todas líneas negras
xlabel('\sl Trayectoria en la zona de Brillouin','FontSize', 10,'Color', 'k');
ylabel('\sl Energía (eV)','FontSize', 10,'Color', 'k');
title('');
set(gca, 'XTick', [])
xlim([1 1000]); % ajustar al rango típico
ylim([-10 10]); % ajustar al rango típico
grid on;