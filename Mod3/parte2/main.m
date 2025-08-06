close all; clear;
pkg load control
addpath("./functions")

% MODELO
[M,K,r] = femPala();
N = length(r); % Cantidad de elementos
Ne = N-1; % Cantidad de elementos
L = diff(r); % Longitudes de los elementos
phi = [1 2 3 4]';
phi_int = [1 1/2 1/3 1/4]';
H_ = [1 0 -3 2; 0 1 -2 1; 0 0 3 -2; 0 0 -1 1];

% CONSTANTES
rhoAir = 1.204; % kg/m^3
CLa = 2*pi;
heli_mass = 450; % kg
w = 550 * (2*pi/60); % Velocidad angular de la pala (rad/s)
c = 0.186; # Cuerda de la pala [m]
g = 9.81; % m/s^2

gamN_b = (1/2)*rhoAir*CLa*c.*(w*r).^2; % Distribución de fuerza 
gamN_yd = (1/2)*rhoAir*CLa*c.*(w*r); % Distribución de fuerza 

phiphi = phi + phi';
phiphi_int = 1 ./(phiphi -1);

F = zeros(2*(Ne+1), 1);
FD = zeros(2*(Ne+1), 2*(Ne+1));
for i = 1:Ne
  le = L(i); % Longitud del elemento
  Dl = [1 0 0 0;
        0 le 0 0;
        0 0 1 0;
        0 0 0 le];
  H = Dl * H_; % Matriz de interpolación del elemento

  gaml_b = (gamN_b(i) + gamN_b(i+1))/2; % Considero distribucion constante como promedio de la Fza en los nodos
  gaml_yd = (gamN_yd(i) + gamN_yd(i+1))/2; % Considero distribucion constante como promedio de la Fza en los nodos

  Fe = gaml_b * le * H * phi_int;
  FeD = gaml_yd * (H * phiphi_int * H')*le; %TODO: CHECK le

  from = 2*i-1;
  to = 2*(i+1);
  F(from:to) = F(from:to) + Fe;
  FD(from:to, from:to) = FD(from:to, from:to) + FeD;
endfor

% CONDICIONES DE BORDE
K(2, :) = [];
K(:, 2) = [];

M(2, :) = [];
M(:, 2) = [];
M(1,1) = M(1,1) + heli_mass/2; % Agregamos la masa del helicóptero

F(2) = [];
FD(2, :) = [];
FD( :,2) = [];

K_ = M\K + eye(size(K))*(w^2);
[V,Lambda] = eig(K_);
Lambda = diag(Lambda);
[~, LIdx] = sort(Lambda);



% Coordenadas modales
F_ = M\F;
P = V\F_;
D_ = M\FD;
D_m = V'*D_*V ;

% Matriz C de salida
C = [-K_(1, :)*V, -D_(1,:)*V];


# función para la variación del paso colectivo

t = [0:0.01:7];
Tau = 0.4;
beta_final = 10*pi/180;
beta = beta_final*(1-e.^(-t/ Tau)); 


figure(1)
plot(t, beta * 180/pi, 'LineWidth', 2);
grid on
ylabel("Cambio de Paso colectivo [grados]");
xlabel("tiempo [s]");
title('Cambio del paso colectivo \Delta\beta(t)');
ylim([0 max(beta * 180/pi)*1.5]);

for nr = 1:4

  % Modelo de orden Reducido

  % Modelo de orden reducido
  idx_a = LIdx(1:nr);
  idx_b = LIdx(nr+1:end);

  L_a = diag(Lambda(idx_a));
  L_b = diag(Lambda(idx_b));
  P_a = P(idx_a, :);
  P_b = P(idx_b, :);
  C_a = C(:, idx_a);
  C_b = C(:, idx_b);
  D_a = D_m(idx_a, idx_a);

  % Modelo de estados
  A_r = [zeros(nr)   eye(nr)
      -L_a         -D_a];

  C_r = [C_a zeros(size(C_a))];

  B_r = [zeros(size(P_a))
      P_a];

  D_r = C_b*inv(L_b)*P_b;

  Mss = ss(A_r, B_r, C_r, D_r);

  [y,t] = lsim(Mss, beta, t);

  figure(2)
  plot(t, 1 + y/g, [";nr = " num2str(nr) " ;"], 'LineWidth', 2);
  grid on
  ylabel("Factor de carga [g]");
  xlabel("tiempo [s]");
  title('Factor de carga g(t)');
  hold on;
end