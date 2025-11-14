 syms W_1L W_R2 epsilon_1 epsilon_2 Delta E q h t real


% Definir la matriz H (Hamiltoniano BdG)
H = [ epsilon_1,     t,      0,     Delta;
           t, epsilon_2,    -Delta,      0;
           0,       -Delta,   -epsilon_1,   -t;
           Delta,     0,        -t,    -epsilon_2 ];

%disp(H);
%gamma= 1; 
gamma= 0.1; 

% Matriz de acoplo W
W = sqrt(gamma/(2*pi))*[
    1,  0, 0, 0; 
    0,  1, 0, 0; 
    0,  0, -1, 0; 
    0,  0, 0, -1];
%disp(W);

% Conjugado transpuesto de W
W_dagger = W';
%disp(W_dagger);

%Matrices identidad y matriz de E. NO ESTOY SEGURO DE SI ES -E O E
I2 = eye(2);
I4 = eye(4);
matrizE = [E, 0, 0, 0;
    0, E,  0, 0;
    0, 0, E, 0;
    0, 0, 0, E];

% B = W*W†
B = W * W_dagger;
A = matrizE - H + 1i*pi*B;
% Inversa de A. NO USAR simplify MEJOR
A_inv = (inv(A));

% Matriz de scattering S. NO USAR simplify MEJOR
S = (I4 - 2*pi*1i * W_dagger * A_inv * W);
%disp(S);

%SELECCIÓN DE DELTA.
%Delta_val = 8*gamma;
Delta_val = 1;

%--------------------------------CÁLCULO DE G------------------------------

%Seleccionamos los elementos S21 y S41: 

tLLhe = S(3,1); %tRL = tLR
tRLee = S(2,1); %aRL = aLR
tRLhe = S(4, 1); %aL

%----------------------------------------
conjtLLhe = conj(tLLhe);
conjtRLee = conj(tRLee);
conjtRLhe = conj(tRLhe);

valor_absoluto_tLLhe = (conjtLLhe)*tLLhe;
valor_absoluto_tRLee = (conjtRLee)*tRLee;
valor_absoluto_tRLhe = (conjtRLhe)*tRLhe;

abs_tLLhe = abs(tLLhe)^2;
abs_tRLee = abs(tRLee)^2;
abs_tRLhe = abs(tRLhe)^2;
%----------------------------------------


%Fórmula de Dourado
G_sym = (q^2/h)*(2*conj(tLLhe)*tLLhe+conj(tRLee)*tRLee+conj(tRLhe)*tRLhe);
disp('G_sym');
disp(G_sym);

title('Cadena de Kitaev: (Método de Dourado)');

%-----------------------------GRÁFICA--2D----------------------------------

%subplot(2,1,2) % Primer cuadro arriba

e1 = linspace(-4.5, 4.5, 100);
e2 = linspace(-4.5, 4.5, 100);
[X1, Y1] = meshgrid(e1, e2);

%Aquí al sustituir el valor de "E" estamos seleccionando el valor de q*V_L
G_subs2 = subs(G_sym, ...
    [ q, h, pi, t, Delta, E], ...
   [1.602e-19, 6.626e-34, 3.141592653589793, Delta_val, Delta_val, 0]);
%disp(G_subs2);

GLR2D = matlabFunction(G_subs2, "File", "GLLmyfilekitaev2acoplosABREU2D", "Vars", [epsilon_1, epsilon_2]);
%disp(GLR2D);

GLR2D_vals = (arrayfun(@(epsilon_1, epsilon_2) GLLmyfilekitaev2acoplosABREU2D(epsilon_1, epsilon_2), X1, Y1))/(((1.602e-19)^2/6.626e-34)*(gamma/Delta_val));
%disp('Valores grafica 2D');
disp((GLR2D_vals));
disp(max(real(GLR2D_vals)));

%Graficamos
%OJO, AUNQUE LA PARTE IMAGINARIA SEA CERO, HAY QUE FORZAR CON real()
surfc(e1/Delta_val, e2/Delta_val, real(GLR2D_vals));
view(2);
shading interp;

%Leyenda
text(0.05, 0.93, '\Gamma = 0.1, \Delta = t = 1, eV_{L} = 0', ...
    'Units', 'normalized', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');
%xticks([-1 -0.5 0 0.5 1]); 
%yticks([-1 -0.5 0 0.5 1]);
xticks([-4 -2 0 2 4]); 
yticks([-4 -2 0 2 4]);
xlim([-4.5 4.5]);  % Ajusta según tus ticks o rango deseado
ylim([-4.5 4.5]);

xlabel('\sl \epsilon_{1} (\Delta)', 'FontSize', 15,'Color', 'k');
ylabel('\sl \epsilon_{2} (\Delta)', 'FontSize', 15,'Color', 'k');
ax2 = gca; % guarda handle

ax2.TickDir = 'out'; % 'in' (por defecto) o 'out'
c = colorbar;
%clim([0, 16]);
colormap(jet);
%c.Ticks = [-7 -5];
%title('Cadena de Kitaev: Dos electrodos (Método de Abreu)');
c.Label.String = 'G_{LL} (\Gamma/\Delta) (e^2/h)';
c.Label.FontAngle = 'italic';
c.Label.FontName = 'Times New Roman';
c.FontSize = 15;

%Ajustar tamaño de gráficas
%set(ax1, 'Position', [0.21 0.6 0.6 0.3]); % [x y width height]
%set(ax2, 'Position', [0.18 0.1 0.6 0.4]);
%set(ax2, 'Position', get(ax2, 'Position'));

% Ajusta colorbar a gráfica
pos_ax2 = ax2.Position;
pos = c.Position;

%Características colorbar
pos(2) = pos_ax2(2);       
pos_cb(4) = pos_ax2(4);      
pos(3) = pos(3) * 0.5;      
grid off;



