clear; clc; close all;
pkg load control

V0=120*1000/3600;
g=9.81;


[t, u, f, St] = getTurbulence(V0);
w = 2*pi*f;

[G, mag] = femRigido(V0, w);


Sy=St.*mag.^2;

y=lsim(G,u,t);
n=(y+g)/g;

figure(1);
plot(f, mag.^2, "LineWidth", 2);
title('Respuesta en frecuencia');
xlabel('Frecuencia [Hz]');
ylabel('|g|^2') ;
grid on;

figure(2);
subplot(2,1,1);
plot(f, Sy, "LineWidth", 2);
title('PSD Salida');
xlabel('Frecuencia [Hz]');
ylabel('Sy');
grid on;

figure(2);
subplot(2,1,2);
plot(t,n, "LineWidth", 2);
grid on;
title("Factor de carga");
ylabel("n [g]");