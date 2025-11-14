syms W_1L W_R2 epsilon_1 epsilon_2 Delta pi E q h t  ...
tLRee tLRhe conjtLRee conjtLRhe ...
valor_absoluto2_tLRee valor_absoluto2_tLRhe G_sym real

%Symbolic Toolbox tiene más precisión que los 16 dígitos por defecto de 
% Matlab y Numpy. Subimos a 25 dígitos de precisión para evitar 
%errores numéricos. 

%Parece que el valor tope de la conductancia cuando el número de dígitos de
%vpa tiende al máximo es del orden de 10^-17. Para vpa de 2^26, sigue 
%dando 10^17. 

%Para valores bajos del dígito de vpa, el tope de la conductancia va
%aumentando, y se pierde precisión. 

%OJO: el número máximo de dígitos de vpa no está especificado en la wiki de
%Matlab, depende del sistema operativo también. En este caso, exceder
%2^26 hará que salte un mensaje por consola

%Este código también funciona, sin embargo, en el caso de no aplicar
%la corrección del infinitésimo, queda GLR = cte. = 0. 

%En teoría, este código es más preciso que GLR_VLneqVR_VPA,
%puesto que forzamos precisión simbólica incluso antes 
%del cálculo de S. 

%Meter vpa en variables
Delta = vpa(1, 25); 
t = vpa(1, 25); 
gamma = vpa(0.1, 25); 
pi = vpa(3.141592653589793, 25); 
E = vpa(0, 25); 
i = vpa(sym(1i), 25);

% Definir la matriz H (Hamiltoniano BdG)
H = [ epsilon_1,     t,      0,     Delta;
           t, epsilon_2,    -Delta,      0;
           0,       -Delta,   -epsilon_1,   -t;
           Delta,     0,        -t,    -epsilon_2 ];


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

%Matrices identidad y matriz de E.
I2 = vpa(sym(eye(2)), 25);
I4 = vpa(sym(eye(4)), 25);

%Sin esta corrección, GLR = cte. = 0
infinitesimo = vpa(sym(1i * 1e-11), 25);
matrizE = (E+infinitesimo)*I4;

%matrizE = [E, 0, 0, 0;
%    0, E,  0, 0;
%    0, 0, E, 0;
%    0, 0, 0, E];

B = W * W_dagger;
A = matrizE - H + i*pi*B;
% Inversa de A. Evitar simplify, puede ser fuente de error
A_inv = inv(A);

% Matriz de scattering S. NO USAR simplify MEJOR
S = (I4 - 2*pi*1i * W_dagger * A_inv * W);
%disp(S);

%--------------------------------CÁLCULO DE G------------------------------

%Seleccionamos los elementos S21 y S41: 

 tLRee = S(1,2);
 tLRhe = S(3,2); 

 %Elementos Dourado
 %tLRee = S(2,1);
 %tLRhe = S(4,1); 

 dif1 = simplify(conj(S(1,2))*S(1,2)-conj(S(2,1))*(S(2,1)));
 dif2 = simplify(conj(S(2,1))*S(2,1)-conj(S(4,1))*(S(4,1)));

 %--------------------------------------------
 conjtLRee = conj(tLRee);
 conjtLRhe = conj(tLRhe);

 valor_absoluto2_tLRee = ((conjtLRee)*tLRee);
 valor_absoluto2_tLRhe = ((conjtLRhe)*tLRhe);

 abs_tLRee = abs(tLRee)^2;
 abs_tLRhe = abs(tLRhe)^2;
%----------------------------------------------

%Fórmula de Dourado
%Forzamos simplify para evitar partes imaginarias

%G_sym = -(valor_absoluto2_tLRee - valor_absoluto2_tLRhe);

G_sym = -(valor_absoluto2_tLRee-valor_absoluto2_tLRhe);

%G_sym = ((-q^2/h)*(conj(tLRee)*(tLRee)-conj(tLRhe)*(tLRhe)));
%disp('G_sym');
%disp(G_sym);

%-----------------------------GRÁFICA--2D----------------------------------

%subplot(2,1,1) % Primer cuadro arriba

e1 = linspace(-4.5, 4.5, 100);
e2 = linspace(-4.5, 4.5, 100);
[X1, Y1] = meshgrid(e1, e2);

GLR2D = matlabFunction(G_sym, "File", "XXxxGLRmyfilekitaev2acoplosABREU2D", "Vars", [epsilon_1, epsilon_2]);
%disp(GLR2D);

GLR2D_vals = vpa(arrayfun(@(epsilon_1, epsilon_2) XXxxGLRmyfilekitaev2acoplosABREU2D(sym(epsilon_1), sym(epsilon_2)), X1, Y1, 'UniformOutput', false), 25);
%disp('Valores grafica 2D');
%disp((GLR2D_vals));
disp(class(GLR2D_vals));

%Convertir de array cell a array
%GLR2D_vals_array_double = double(reshape([GLR2D_vals{:}], size(GLR2D_vals)));

%Aqui es un array simbolico, no un array cell, luego usamos double
%directamente

GLR2D_vals_array_double = double(GLR2D_vals);

%Pasamos a escala Gamma/Delta
GLR2D_vals1 = GLR2D_vals_array_double/double(gamma/Delta);

%Graficamos
%OJO, AUNQUE LA PARTE IMAGINARIA SEA CERO, HAY QUE FORZAR CON real()
surfc(e1, e2, real((GLR2D_vals1)));
view(2);
shading interp;


%max_abs_val = max(((GLR2D_vals(:))));
%min_abs_val = min((GLR2D_vals(:)));
%disp(max_abs_val);
%disp(min_abs_val);

%Leyenda
%text(0.05, 0.93, '\Gamma = 1, \Delta = t = 8\Gamma, eV_{R} = 0.1\rm ', ...
%    'Units', 'normalized', ...
%    'FontSize', 10, ...
%    'BackgroundColor', 'w', ...
%    'EdgeColor', 'k');
%text(0.05, 0.93, '\Gamma = 1, \Delta = 8\Gamma, t = 4\Gamma, eV_{R} = 0.0001\rm ', ...
%    'Units', 'normalized', ...
%    'FontSize', 10, ...
%    'BackgroundColor', 'w', ...
%    'EdgeColor', 'k');
%text(0.05, 0.93, '\Gamma = 1, \Delta = t = 8\Gamma, eV_{R} = 0\rm ', ...
%    'Units', 'normalized', ...
%    'FontSize', 10, ...
%    'BackgroundColor', 'w', ...
%    'EdgeColor', 'k');
text(0.05, 0.93, '\Gamma = 0.1, \Delta = t = 1, eV_{R} = 0\rm ', ...
    'Units', 'normalized', ...
    'FontSize', 10, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k');
%xticks([-1.5 -1 -0.5 0 0.5 1 1.5]); 
%yticks([-1.5 -1 -0.5 0 0.5 1 1.5]);
xticks([-4 -2 0 2 4]); 
yticks([-4 -2 0 2 4]);
%xlim([-1.5 1.5]);  % Ajusta según tus ticks o rango deseado
%ylim([-1.5 1.5]);
xlim([-4.5 4.5]);  % Ajusta según tus ticks o rango deseado
ylim([-4.5 4.5]);

%LÍMITES BARRA DE COLOR
%clim([-0.000000000000002, 0.000000000000002]);

xlabel('\sl \epsilon_{1} (\Delta)', 'FontSize', 15,'Color', 'k');
ylabel('\sl \epsilon_{2} (\Delta)', 'FontSize', 15,'Color', 'k');
ax2 = gca; % guarda handle

ax2.TickDir = 'out'; % 'in' (por defecto) o 'out'
c = colorbar;
colormap(jet);
c.Label.String = 'G_{LR} (\Gamma/\Delta) (e^2/h)';
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

% Opcional: mover a la derecha del eje
%pos(1) = pos_ax2(1) + pos_ax2(3) + 0.02;
%c.Position = pos;
grid off;



