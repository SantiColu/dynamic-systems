close all; clear;
pkg load control
pkg load signal

load("aceleracion_g.mat");
accel = y;
load("fem_lng.mat");

N = length(X);
K = K(2:end, 2:end);
M = M(2:end, 2:end);

K_ = M\K;
F = K_(:,1);
C = zeros(1, N-1);
C(end) = 1;

K_ = M\K;
F_ = M\F;
[V,L] = eig(K_, 'vector');
[~, LIdx] = sort(L);

%Coordenadas modales
C = C*V;
P = V\F_;

[V,L] = eig(K_, 'vector');
[~, LIdx] = sort(L);

%Coordenadas modales
C = C*V;
P = V\F_;

% Modelo de orden Reducido
nr = 5;

idx_a = LIdx(1:nr);
idx_b = LIdx(1:nr);

L_a = diag(L(idx_a));
L_b = diag(L(idx_b));
P_a = P(idx_a, :);
P_b = P(idx_b, :);
C_a = C(:, idx_a);
C_b = C(:, idx_b);


% Modelo de estados
A = [zeros(nr)   eye(nr)
    -L_a         zeros(nr)];

B = [zeros(size(P_a))
    P_a];

C = [C_a zeros(size(C_a))];

D = C_b*inv(L_b)*P_b;

Mss = ss(A, B, C, D);
G = tf(Mss);


% ---------------------------

fs = 5000;
ts = 1/fs;
t = 0:ts:(length(accel)-1)*ts;
u = accel;


output = lsim(G, u, t);



nw = 2^6;
Su = pwelch(u, nw, 0.5, nw, fs);
Sy = pwelch(output, nw, 0.5, nw, fs);
Df = fs/nw;
f = (1:length(Su))*Df;

tM = length(accel)/5;
figure(1)
subplot(2,1,1)
plot(t(1:tM), u(1:tM));
ylabel(" Aceleración [m/s²]")
xlabel('Tiempo [s]');
title("Entrada");
grid on;

subplot(2,1,2)
plot(t(1:tM), output(1:tM));
ylabel(" Aceleración [m/s²]")
xlabel('Tiempo [s]');
title("Salida");
grid on;



figure(2)
subplot(1,2,1)
semilogx(f, Su, "LineWidth", 2);
ylabel('PSD');
xlabel('Frecuencia [hz]');
title("Entrada");
grid on;

subplot(1,2,2)
semilogx(f, Sy, "LineWidth", 2);
ylabel('PSD');
xlabel('Frecuencia [Hz]');
title("Salida");
grid on;


% Comparacion diferentes tamaño de bloques
for i = 5:8 
  nw = 2^i;
  Su = pwelch(u, nw, 0.5, nw, fs);
  Sy = pwelch(output, nw, 0.5, nw, fs);
  Df = fs/nw;
  f = (1:length(Su))*Df;

  figure(3)
  subplot(1,2,1)
  semilogx(f, Su, ["; Bloque 2^", int2str(i), ";"], "LineWidth", 2);
  ylabel('PSD');
  xlabel('Frecuencia [hz]');
  title("Entrada");
  grid on;
  hold on;

  subplot(1,2,2)
  semilogx(f, Sy, ["; Bloque 2^", int2str(i), ";"], "LineWidth", 2);
  ylabel('PSD');
  xlabel('Frecuencia [Hz]');
  title("Salida");
  grid on;
  hold on;
endfor
