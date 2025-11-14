
Delta_MP = 1;
mu_MP = linspace(-4, 4, 101);
tau_MP = linspace(0, 4, 101);

%--------------------------------------------------------------------------

%GRÁFICA SUPERIOR
[X1, Y1] = meshgrid(tau_MP, mu_MP);
M = (1./sqrt(1+Y1.^2));

subplot(2,1,1) % Primer cuadro arriba
surfc(tau_MP/Delta_MP, mu_MP/Delta_MP, M);
view(2);
shading interp;
xticks([0 1 2 3 4]);
yticks([-4 -3 -2 -1 0 1 2 3 4]);
zticks([0 1]);
xlabel('\sl t (\Delta)', 'FontSize', 20,'Color', 'k');
ylabel('\sl \mu (\Delta)', 'FontSize', 20,'Color', 'k');
c = colorbar;
clim([0 1]);
c.Ticks = [0 1];
c.Label.String = '|M|';
c.Label.FontAngle = 'italic';
c.Label.FontName = 'Times New Roman';
c.FontSize = 20;
%----------------------Customizamos la escala de color---------------------
% Número de colores
n = 256;

% Primer subplot: 0 → azul oscuro, 1 → amarillo
ax1 = subplot(2,1,1);

% Fijar escala de color de 0 a 1
caxis([0 1]);

% Colormap: azul oscuro → amarillo
blue = [0, 0, 1];
yellow = [1, 1, 0.2];
cmap1 = [linspace(blue(1), yellow(1), n)', ...
         linspace(blue(2), yellow(2), n)', ...
         linspace(blue(3), yellow(3), n)'];
colormap(ax1, cmap1);

ax1.TickDir = 'out'; % 'in' (por defecto) o 'out'

ax1.FontSize = 18;         % cambia tamaño de los números de X e Y

%-------------------------------------------------------------------------
box on;
grid on;
xLimits = xlim;
yLimits = ylim;
% Dibuja el marco negro
%rectangle('Position', [xLimits(1), yLimits(1), ...
%                       diff(xLimits), diff(yLimits)], ...
%
%           'EdgeColor', 'k', 'LineWidth', 2);

% title('$\mathrm{MP}$', ...
%       'Interpreter', 'latex', ...
%       'FontSize', 15);
hold on;

%--------------------------------------------------------------------------

%GRÁFICA INFERIOR
[X2, Y2] = meshgrid(tau_MP, mu_MP);
deltaE = sqrt(Delta_MP.^2+Y2.^2)-X2;
subplot(2,1,2) % Segundo cuadro abajo
surfc(tau_MP/Delta_MP, mu_MP/Delta_MP, deltaE/Delta_MP);
view(2);
shading interp;

xticks([0 1 2 3 4]);
yticks([-4 -3 -2 -1 0 1 2 3 4]);
zticks([-4 -2 2 4]);
xlabel('\sl t (\Delta)', 'FontSize', 20,'Color', 'k');
ylabel('\sl \mu (\Delta)', 'FontSize', 20,'Color', 'k');
c = colorbar;
clim([-4 4]);
c.Ticks = [-4 -2 2 4];
c.Label.String = '\deltaE_0 (\Delta)';
c.Label.FontAngle = 'italic';
c.Label.FontName = 'Times New Roman';
%c.FontWeight = 'bold';
c.FontSize = 20;

%----------------------Customizamos la escala de color---------------------
ax2 = subplot(2,1,2);

% Fijar escala de color de -4 a 4
caxis([-4 4]);

% Colormap: púrpura → blanco → rojo
purple = [0.5, 0, 0.5];
white = [1, 1, 1];
red = [1, 0, 0];

half1 = [linspace(purple(1), white(1), n/2)', ...
         linspace(purple(2), white(2), n/2)', ...
         linspace(purple(3), white(3), n/2)'];

half2 = [linspace(white(1), red(1), n/2)', ...
         linspace(white(2), red(2), n/2)', ...
         linspace(white(3), red(3), n/2)'];

cmap2 = [half1; half2];
colormap(ax2, cmap2);
ax2.FontSize = 18;         % cambia tamaño de los números de X e Y
ax2.TickDir = 'out'; % 'in' (por defecto) o 'out'
%-------------------------------------------------------------------------

box on;
grid on;

xLimits = xlim;
yLimits = ylim;

% Dibuja el marco negro
%rectangle('Position', [xLimits(1), yLimits(1), ...
%                      diff(xLimits), diff(yLimits)], ...
%         'EdgeColor', 'k', 'LineWidth', 2);
% title('$\mathrm{\delta E_{0}} = E_{0}^{odd}-E_{0}^{even}$','Interpreter', 'latex', 'FontSize', 15);
hold off;


