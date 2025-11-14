clear all;
close all;
clc;

% Parámetros
a = 2.46;
Es = -8.37;
Epx = 0;
Epy = 0;
Epz = 1;   
Vsssigma = -5.73;
Vspsigma = 6.05;
Vpppi = -3.07;
Vppsigma = 5.62;


%Definimos la matriz por bloques

        H = zeros(60, 60);
        for i = 1:15

            h = i*4-3;
            g = i*4-3;

            for k = 1:3
            H(g+k, h) = Vpppi;
            end
        end 


        for i = 1:15

            h = i*4-3;
            g = i*4-3;

            for k = 1:3
            H(g, h+k) = Vpppi;
            end
        end 


        for i = 1:60
            for j = 1:60
                if i == j
                   H(i, j) = Epz;
                end
            end  
        end
           
        disp(H);
        
        % Calcular autovalores
        E = real(eig(H));

        %Mostrar autovalores
        disp(E);
    
%Visualizar el espectro energético

x1 = -0.5;
x2 = 0.5;

figure;
hold on;
for i = 1:length(E)
    plot([x1 x2], [E(i) E(i)], 'b-', 'LineWidth', 2);
end

set(gca, 'xtick', []);
ylabel('\sl Energía (eV)','FontSize', 10,'Color', 'k');
title('Niveles energéticos del C60');
xlim([-1 1]);
ylim([min(E)-0.5, max(E)+0.5]);
