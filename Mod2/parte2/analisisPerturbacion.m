close all; clear;
addpath("./functions")

[M,K,r,c] = femPala();

w = 88 * (2*pi/60);         % Velocidad angular de la pala (rad/s)
T = 2*pi/w;
[t, v] = modelWind(30, 11.3, max(r), w);

% Plot de la velocidad del viento
dt = t(2) - t(1); % Paso de tiempo
figure(1);
plot(t(1:3*T/dt), v(1:3*T/dt), 'b', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('Velocidad del viento (por variacion vertical) [m/s]');
grid on;

% Aproximación de la velocidad del viento
[v_aprox, m] = fftWind(t, v, true);

% Comparar
figure
plot(t(1:3*T/dt), v_aprox(1:3*T/dt), 'r--', 'LineWidth', 1.5);
hold on;
plot(t(1:3*T/dt), v(1:3*T/dt), 'b', 'LineWidth', 1.5);
legend('Aprox armónica dominante', 'Original');
xlabel('Tiempo [s]'); ylabel('V(t)');
title('Comparación entre señal original y su armónico dominante');
grid on;
