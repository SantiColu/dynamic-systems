close all; clear;

addpath("./functions")


[M,K] = femGeneral();

% Análisis modal
K_ = M\K;
[V,D] = eig(K_, 'vector');
[D_s, dIdx] = sort(D);
f = real(sqrt(D_s)/(2*pi)); %frecuencias de los modos en hz

% Nos quedamos solo con algunas
f = f(1:3*4);

% Mostrar frecuencias
disp('Frecuencias naturales [Hz]:');
disp(f);

% Mostrar frecuencias
disp('Frecuencias naturales [RPM]:');
disp(f*60);

% Graficar frecuencias
plotProhibidas(1, f, 'Frecuencia (Hz)', 'Zonas prohibidas', 'Zonas de frecuencia prohibidas unidas (±10%)')
plotProhibidas(2, f *60, 'Frecuencia (RPM)', 'Zonas prohibidas', 'Zonas de frecuencia prohibidas unidas (±10%)')



