close all; clear;
pkg load control
addpath("./functions")

% CONSTANTES
rhoAir = 1.204; % kg/m^3
CLa = 2*pi;
alfa_n = 6.5*pi/180;

% MODELO
[M,K,r,c] = femPala();
N = length(r); % Cantidad de elementos
Ne = N-1; % Cantidad de elementos
L = diff(r); % Longitudes de los elementos
phi_int = [1 1/2 1/3 1/4]';
H_ = [1 0 -3 2; 0 1 -2 1; 0 0 3 -2; 0 0 -1 1];

%Condicion nominal
wn =  65 * (2*pi/60);       % Velocidad nominal angular de la pala (rad/s)
Vn = 11.3;                % Velocidad nominal del viento (m/s)
phi_n = atan(wn*r/Vn);
beta = alfa_n + phi_n;


% Matriz C de salida
C = zeros(1, 2*N);
C(3) = K(3,2);
C(4) = K(4,2);

% Empotramos la raiz de la pala
M = M(3:end, 3:end);
K = K(3:end, 3:end);
C = C(3:end);

K_ = M\K;
[Vect,Lambda] = eig(K_, 'vector');
[~, LIdx] = sort(Lambda);

% w_rng = logspace(0, 4, 1000); % Rango de frecuencias (rad/s)
w_rng = linspace(0, 10, 100); % Rango de frecuencias (rad/s)
M0_w = zeros(length(w_rng), 1);
M_w = zeros(length(w_rng), 1);
g_w = zeros(length(w_rng), 1);
d_w = zeros(length(w_rng), 1);
for k = 1:length(w_rng)
  w = w_rng(k); %* (2*pi/60);         % Velocidad angular de la pala (rad/s)

  % Perturbación del viento
  [t, v] = modelWind(30, 11.3, max(r), w);
  [_, m, a0] = fftWind(t, v);

  phi = atan(w*r/Vn);
  alfa = beta - phi;
  gamN = (1/2)*rhoAir*CLa*c.*alfa.*sin(phi); % Distribución de fuerza %TODO: GENERAR LA DISTRIBUCION DE FZA

  F = zeros(2*(Ne+1), 1);
  Fs = zeros(2*(Ne+1), 1);
  for i = 1:Ne
    le = L(i); % Longitud del elemento
    Dl = [1 0 0 0;
          0 le 0 0;
          0 0 1 0;
          0 0 0 le];
    H = Dl * H_; % Matriz de interpolación del elemento

    gaml = (gamN(i) + gamN(i+1))/2; % Considero distribucion constante como promedio de la Fza en los nodos
    re = (r(i) + r(i+1))/2; % Radio medio del elemento

    Fe = gaml * le * H * phi_int;

    from = 2*i-1;
    to = 2*(i+1);
    F(from:to) = F(from:to) + Fe;
    Fs(from:to) = Fs(from:to) + Fe.*(w*re)^2;
  endfor

  % Empotramos la raiz de la pala
  F = F(3:end);
  Fs = Fs(3:end);

  % Momento estatico
  q0 = K\Fs;
  M0 = C*q0;

  M0_w(k) = M0; % Momento en la raiz de la pala

  % Coordenadas modales
  F_ = M\F;
  P = Vect\F_;
  C_ = C*Vect;

  % Modelo de orden Reducido
  nr = 5;

  % Modelo de orden reducido
  idx_a = LIdx(1:nr);
  idx_b = LIdx(1:nr);

  L_a = diag(Lambda(idx_a));
  L_b = diag(Lambda(idx_b));
  P_a = P(idx_a, :);
  P_b = P(idx_b, :);
  C_a = C_(:, idx_a);
  C_b = C_(:, idx_b);

  chi = 0.05;

  % Modelo de estados
  A_r = [zeros(nr)   eye(nr)
      -L_a         -2*chi*sqrt(L_a)];

  C_r = [C_a zeros(size(C_a))];

  B_r = [zeros(size(P_a))
      P_a];

  D_r = C_b*inv(L_b)*P_b;

  Mss = ss(A_r, B_r, C_r, D_r);

  % Analisis en frecuencia
  G = tf(Mss);
  [g0, ~] = bode(G, 0);
  [g, d] = bode(G, w);
  M_max =  -g0*a0 + g*m + M0;
  M_w(k) = M_max;
  g_w(k) = g;
  d_w(k) = d;
endfor


figure(1);
plot(w_rng * 30/pi, M0_w, 'b', 'LineWidth', 1.5);
xlabel('Frecuencia [RPM]');
ylabel('Momento estatico [Nm]');
grid on;

figure(2);
plot(w_rng * 30/pi, M_w, 'b', 'LineWidth', 1.5);
xlabel('Frecuencia [RPM]');
ylabel('Momento total [Nm]');
grid on;

figure(3);
subplot(211);
semilogx(w_rng, 20*log(g_w), 'b', 'LineWidth', 1.5)
xlabel('Frecuencia [rad/s]');
ylabel('Ganancia [db]');
title('Ganancia');
grid on;

subplot(212);
semilogx(w_rng, d_w, 'b', 'LineWidth', 1.5)
xlabel('Frecuencia [rad/s]');
ylabel('Desfasaje [grados]');
title('Desfasaje');
grid on;