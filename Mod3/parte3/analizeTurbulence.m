clear; close all;
pkg load control
addpath("./functions")


V0 = 120/3.6; % Velocidad de vuelo [m/s]

[t, u, f, St] = getTurbulence(V0);

figure(1);
subplot(2,1,1);
plot(f, St, "LineWidth", 2);
title('PSD Turbulencia');
ylabel('St');
grid on;

subplot(2,1,2);
plot(t,u, "LineWidth", 2);
grid on;
title("Turbulencia atmosférica");
xlabel("t [s]");
ylabel("w_w [m/s]");