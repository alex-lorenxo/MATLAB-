syms W_1L W_R2 epsilon_1 epsilon_2 Delta E q h t real

% Definir la matriz H (Hamiltoniano BdG)
H = [ epsilon_1,     t,         0,     Delta;
           t, epsilon_2,    -Delta,      0;
           0,   -Delta,   -epsilon_1,   -t;
        Delta,     0,        -t,    -epsilon_2 ];
%disp(H);

%Valor del parámetro de acoplo en wide-band
gamma= 1; 

% Matriz de acoplo W
W = sqrt(gamma/(2*pi))*[1,  0,   0,   0;
     0,   0, -1,   0];
%disp(W);

% Conjugado transpuesto de W
W_dagger = W';
%disp(W_dagger);

%Matrices identidad y matriz de E
I2 = eye(2);
I4 = eye(4);
matrizE = [E, 0, 0, 0;
    0, E,  0, 0;
    0, 0, E, 0;
    0, 0, 0, E];
%La fórmula de Flensberg no lleva un - en el bloque de huecos de la matriz
% E, a diferencia del formalismo de Maiani y 


% B = W*W†
B = W_dagger*W;
A = H - matrizE + 1i*pi*B;
% Inversa de A. NO USAR simplify MEJOR
A_inv = (inv(A));

% Matriz de scattering S. NO USAR simplify MEJOR
S = (I2 + 2*1i*pi*W*A_inv*W_dagger);
%disp(S);

% Submatrices r11, r22, t12, t21
r11 = S(1:1, 1:1);
t21 = S(2:2, 1:1);
t12 = S(1:1, 2:2);
r22 = S(2:2, 2:2);

%disp(r11);
%disp(t21);
%disp(t12);
%disp(r22);

%Valores de la abscisa en eV
E_vals = (1.602e-19)*linspace(-3.1210986e+20, 3.1210986e+20, 1000); 

%--------------------------DIAGONALIZACIÓN DE H---------------------------
%---------------------DIAGONALIZACIÓN ANALÍTICA---------------------------

[V, D] = eig(H);
%Mostrar autovalores por consola 
%disp(simplify(D));

%-----------------------DIAGONALIZACIÓN NUMÉRICA---------------------------

%Sustituimos valores en el sweet spot  
H_subs = subs(H, [epsilon_1, epsilon_2, t, Delta], [0, 0, 8*gamma, 8*gamma]);
[V_subs, D_subs] = eig(H_subs); 
%Mostrar autovectores y autovalores por consola
%disp(V_subs); 
%disp(D_subs);
%gamma = 1, luego t = Delta = 8. Y deberiamos obtener por consola los
%autovalores 0, 0, 16, -16

%--------------------------------CÁLCULO DE G------------------------------
%Flensberg solo tiene en cuenta la reflexión de Andreev de electrones a, 
%correspondiente al elemento de la matriz de scattering t21
G_sym = (2*q^2/h)*(conj(t12)*t12);

%disp(G_sym);

%Sustitución de parámetros
G_subs = subs(G_sym, ...
    [ q, h, epsilon_1, epsilon_2, pi, t], ...
    [1.602e-19, 6.626e-34, 8*gamma, 0, 3.141592653589793, 8*gamma]);

G = matlabFunction(G_subs, "File", "myfilekitaevFLENSBERG", "Vars", [E, Delta]);

%Selección de Delta
Delta_val = 8*gamma;

% Evaluar la función. Usando arrayfun, generamos el vector G_vals 
%G_vals = arrayfun(@(E) myfilekitaevFLENSBERG(E, Delta_val) + myfilekitaevFLENSBERG(-E, Delta_val), E_vals);
G_vals = arrayfun(@(E) myfilekitaevFLENSBERG(E, Delta_val), E_vals);


% Graficamos la conductancia dividida por el factor (2e^2/h)
plot(E_vals/gamma, (G_vals)/(2*(1.602e-19)^2/6.626e-34), 'LineWidth', 4);
xlabel('\sl eV (\Gamma)', 'FontSize', 50,'Color', 'k');
xlim([-20, 20]);
xticks([-20 -10 0 10 20]);
yticks([0 0.5 1]);
ylabel('\sl G(2e^2/h)', 'FontSize', 50,'Color', 'k');

lgd = legend("\epsilon_{1} = 0, \epsilon_{2} = 0, \Delta = 8\Gamma, t = 8\Gamma", ...
    "\epsilon_{1} = 4\Gamma, \epsilon_{2} = 0, \Delta = 8\Gamma, t = 8\Gamma", ...
    "\epsilon_{1} = 8\Gamma, \epsilon_{2} = 0, \Delta = 8\Gamma, t = 8\Gamma"); 
lgd.FontSize = 30; 


lgd = legend('Location','eastoutside','Box','on'); 
lgd.AutoUpdate = 'off';  

ax = gca;          % obtiene el eje actual
ax.FontSize = 30;  % cambia el tamaño de los números de los ticks


title('(a)');
hold on; 
grid on;
