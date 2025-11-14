% Leer los datos del archivo table.txt sin modificar el nombre del archivo
table = readtable('table.txt', 'FileType', 'text', 'VariableNamingRule', 'preserve');

% Extraer la columna de tiempo
time = table{:, 1}; % Suponemos que la primera columna es '# t (s)'

% Número de regiones por rama
num_regions_per_branch = 45; % Ahora hay 45 regiones por Entrada 1, Entrada 2 y Salida

% Crear un vector de regiones en nanómetros
regions_input = (0:num_regions_per_branch - 1) * 8; % Distancia en nm entre los centros
regions_output = regions_input;

% Inicializar matrices para Entrada 1, Entrada 2 y Salida
mx1 = zeros(length(time), num_regions_per_branch);
my1 = zeros(length(time), num_regions_per_branch);

mx2 = zeros(length(time), num_regions_per_branch);
my2 = zeros(length(time), num_regions_per_branch);

mx3 = zeros(length(time), num_regions_per_branch);
my3 = zeros(length(time), num_regions_per_branch);

% Desplazamiento inicial debido a columnas adicionales
base_col_offset = 5;

% Extraer las componentes de la magnetización
% Entrada 1
for i = 1:num_regions_per_branch
    col_idx_x = base_col_offset + 3 * (i - 1); % Índice para la componente x
    col_idx_y = col_idx_x + 1;   % Índice para la componente y
    mx1(:, i) = table{:, col_idx_x};
    my1(:, i) = table{:, col_idx_y};
end


% Entrada 2
for i = 1:num_regions_per_branch
    col_idx_x = base_col_offset + num_regions_per_branch * 3 + 3 * (i - 1);
    col_idx_y = col_idx_x + 1;
    mx2(:, i) = table{:, col_idx_x};
    my2(:, i) = table{:, col_idx_y};
end


% Salida (45 regiones)
for i = 1:num_regions_per_branch
    col_idx_x = base_col_offset + 2 * num_regions_per_branch * 3 + 3 * (i - 1);
    col_idx_y = col_idx_x + 1;
    mx3(:, i) = table{:, col_idx_x};
    my3(:, i) = table{:, col_idx_y};
end


%Definimos el vídeo
obj = VideoWriter('animacion_ondas_de_spin.avi');
obj.Quality = 70;
obj.FrameRate = 60;
open(obj); 

%fig = figure('Position', [200 200 600 600], 'Visible', 'on');

fig = figure('WindowState', 'maximized'); % Maximiza la figura a pantalla completa


% Inicializar el título global con el primer instante de tiempo
t = sgtitle(sprintf('Instante de tiempo: %.2f ps', 0 * 2.5));

for cnt=1:length(time)

% Configurar subgráficos
legend('Interpreter', 'latex'); % Activa LaTeX en la leyenda
% Control
pos1 = [0.12, 0.68, 0.8, 0.2]
subplot('position', pos1);

plot(regions_input, mx1(cnt, :), 'r', 'DisplayName', 'm_x'); hold on;
plot(regions_input, my1(cnt, :), 'g', 'DisplayName', 'm_y');
hold off;
title('Control', 'FontSize', 20);
ylabel({'Magnetización ', 'normalizada'}, 'FontSize', 18);
legend;
grid on;
xlim([min(regions_input), max(regions_input)]);
set(gca, 'XTickLabel', [], 'FontSize', 18) % Elimina los números del eje x pero mantiene las marcas
ylim([-0.2, 0.2]);
yticks(-0.2:0.1:0.2); % Especifica las marcas del eje y cada 0.05 unidades

% Entrada 
pos2 = [0.12, 0.43, 0.8, 0.2]
subplot('Position', pos2);

plot(regions_input, mx2(cnt, :), 'r', 'DisplayName', 'm_x'); hold on;
plot(regions_input, my2(cnt, :), 'g', 'DisplayName', 'm_y');
hold off;
title('Entrada', 'FontSize', 20);
ylabel({'Magnetización ', 'normalizada'}, 'FontSize', 18);
legend;
grid on;
xlim([min(regions_input), max(regions_input)]);
set(gca, 'XTickLabel', [], 'FontSize', 18); % Elimina los números del eje x pero mantiene las marcas
ylim([-0.2, 0.2]);
yticks(-0.2:0.1:0.2); % Especifica las marcas del eje y cada 0.05 unidades

% Salida (45 regiones)
pos3 = [0.12, 0.18, 0.8, 0.2]
subplot('Position', pos3);

plot(regions_output, mx3(cnt, :), 'r', 'DisplayName', 'm_x'); hold on;
plot(regions_output, my3(cnt, :), 'g', 'DisplayName', 'm_y');
hold off;
title('Salida', 'FontSize', 20);
xlabel('Posición (nm)', 'FontSize', 18);
ylabel({'Magnetización ', 'normalizada'}, 'FontSize', 18);
legend;
grid on;
xlim([min(regions_input), max(regions_input)]);
set(gca, 'FontSize', 18);
ylim([-0.2, 0.2]);
yticks(-0.2:0.1:0.2); % Especifica las marcas del eje y cada 0.05 unidades


hold off; 





% Actualizar el título con el instante de tiempo correcto en picosegundos
t.String = sprintf('Instante de tiempo: %.2f ps', (cnt-1) * 2.5);

f = getframe(fig); % Captura la figura oculta
writeVideo(obj,f);


end
close(obj);


