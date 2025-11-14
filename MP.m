 
delta = 1;
mu = linspace(-4, 4, 101);
tau = linspace(0, 4, 101);

[X, Y] = meshgrid(tau, mu);
M = (1./sqrt(1+Y.^2));

surfc(tau/delta, mu/delta, M);
shading interp;
xticks([0 2 4]); 
yticks([-4 -2 0 2 4]);
zticks([0 1]);
%zticklabels({'z = 0', 'z = 1'});
xlabel('\sl \tau (\delta)', 'FontSize', 15,'Color', 'k');
ylabel('\sl \mu (\delta)', 'FontSize', 15,'Color', 'k');
c = colorbar;
clim([0 1]);
c.Ticks = [0 1];
c.Label.String = '|M|';
d = listfonts;
disp(d);
c.Label.FontAngle = 'italic';
c.Label.FontName = 'Times New Roman';
%c.FontWeight = 'bold';
c.FontSize = 15;
box on;
grid off;

xLimits = xlim;
yLimits = ylim;

% Dibuja el marco negro
hold on
rectangle('Position', [xLimits(1), yLimits(1), ...
                       diff(xLimits), diff(yLimits)], ...
          'EdgeColor', 'k', 'LineWidth', 2)

title('$\mathrm{MP}(\mu_1 = \mu_2 = \mu)$', ...
      'Interpreter', 'latex', ...
      'FontSize', 15);
%title('MP (\mu_1 = \mu_2 = \mu)', 'FontSize', 13);











