clear; clc; close all;
pkg load control
addpath("functions");

V0=120*1000/3600;
g=9.81;


[t, u, f, St] = getTurbulence(V0);
w = 2*pi*f;

[G_sl, mag_sl] = femSinLastre(V0, w);
Sy_sl=St.*mag_sl.^2;
n_sl=(lsim(G_sl,u,t)+g)/g;

[G_r, mag_r] = femRigido(V0, w);
Sy_r=St.*mag_r.^2;
n_r=(lsim(G_r,u,t)+g)/g;

[G_cl, mag_cl] = femConLastre(V0, w);
Sy_cl=St.*mag_cl.^2;
n_cl=(lsim(G_cl,u,t)+g)/g;


figure(1);
title('Respuesta en frecuencia');
ylabel('|g|^2') ;
xlabel('Frecuencia [Hz]');
grid on;
hold on;
plot(f, mag_sl.^2, ";Sin Lastre;", "LineWidth", 2);
plot(f, mag_r.^2, ";Cuerpo Rigido;", "LineWidth", 2);
plot(f, mag_cl.^2, ";Con Lastre;", "LineWidth", 2);

figure(2);
subplot(2,1,1);
title('PSD Salida');
xlabel('Frecuencia [Hz]');
ylabel('Sy');
grid on;
hold on;
plot(f, Sy_sl, ";Sin Lastre;", "LineWidth", 2);
plot(f, Sy_r, ";Cuerpo Rigido;", "LineWidth", 2);
plot(f, Sy_cl, ";Con Lastre;", "LineWidth", 2);

figure(2);
subplot(2,1,2);
grid on;
title("Factor de carga");
ylabel("n [g]");
hold on;
plot(t,n_sl, ";Sin Lastre;", "LineWidth", 2);
plot(t,n_r, ";Cuerpo Rigido;", "LineWidth", 2);
plot(t,n_cl, ";Con Lastre;", "LineWidth", 2);
