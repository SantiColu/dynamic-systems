close all; clear;
h0 = 30;
V0 = 11.3;

h = linspace(0, h0 + 2*7.4900, 500);
U = V0 * (h/h0).^0.4;

h1 = h0 - 7.4900;
h2 = h0 + 7.4900;

idx = h >= h1 & h <= h2;
h_fill = h(idx);
U_fill = U(idx);


hold on;


fill([zeros(size(U)), fliplr(U)], [h, fliplr(h)], ...
     [124 228  255]/255, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

pl(1) = fill([zeros(size(U_fill)), fliplr(U_fill)], [h_fill, fliplr(h_fill)], ...
     [255 189 89]/255, 'FaceAlpha', 0.8, 'EdgeColor', 'none',  "DisplayName", "Velocidades de interes");  


xlim()

pl(2) = plot(U, h, 'b', 'LineWidth', 1.5, "DisplayName", "Perfil de Velocidades");

axis equal;

pl(3) = plot(xlim(), [h1 h1], '--', 'color', "#ffbd59", 'DisplayName', 'Altura minima de pala');
pl(4) = plot(xlim(), [h2 h2], '--', 'color', "#ffbd59", 'DisplayName', 'Altura maxima de pala');

xlabel('Velocidad del viento [m/s]');
ylabel('Altura [m]');
title('Perfil de velocidad del viento');
legend(pl)

grid on;
