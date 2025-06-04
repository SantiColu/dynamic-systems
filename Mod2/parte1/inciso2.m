clear; close all;
pkg load control;
addpath("./functions")

% condición de vuelo
h   = 1000; % m
V   = 140 ; % m/s

addpath('../../Mod1/mdl');
mdl = A320_build_model(1, h, V/1.852*3.6, 0, false, false, 0, 0);
g = mdl.g;
rho = mdl.rho;
U = mdl.speed * cosd(mdl.alfa);
W = mdl.mass;
q = mdl.q;
alfa = mdl.alfa*pi/180;

addpath('../../Mod1/wing');
N = 100;
[M,K,x,elements] = modelwing(N,0);

C = zeros(1, 2*N);
C(3) = K(3,2);
C(4) = K(4,2);

% Matriz F
At = sum(elements.area);
CL = W/2*g/(q*At*alfa);

F = zeros(2*N, 1);
for i = 1:N
  F(2*i-1) = q*CL*elements.area(i);
end

% Empotramos la raiz del ala
M = M(3:end, 3:end);
K = K(3:end, 3:end);
C = C(3:end);
F = F(3:end);

% Momento estatico
Fe = F.*alfa;
q0 = K\Fe;
M0 = C*q0;

K_ = M\K;
F_ = M\F;
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

chi = 0.05;

% Modelo de estados
A = [zeros(nr)   eye(nr)
    -L_a         -2*chi*sqrt(L_a)];

B = [zeros(size(P_a))
    P_a];

C = [C_a zeros(size(C_a))];

D = C_b*inv(L_b)*P_b;

stateName = {};
for i = 1:nr
    stateName(end+1) = sprintf('u%d', i);
    stateName(end+1) = sprintf('t%d ', i);
endfor
inputName = {'ww'};
outputName = {"M_raiz"};

Mss = ss(A, B, C, D, 'StateName', stateName, 'InputName', inputName, 'OutputName', outputName);

% Rafaga (1-cos)
H = 9.144;
Uw = 10;
[t,ww] = modelBurst("1-cos", H, Uw, U, 30);
u = ww/(mdl.speed*cos(alfa));
[m, ~] = lsim(Mss, u, t);

figure(1);
subplot(2, 1, 1);
plot(t, u, "LineWidth", 2);
title('Rafaga (1-cos)');
xlabel('t [s]');
ylabel('u [m/s]');
grid on;
hold on;

subplot(2, 1, 2);
plot(t, m + M0, "LineWidth", 2);
xlabel('t [s]');
ylabel('M [Nm]');
title('Respuesta (M_{raiz})');
grid on;
hold on;


[t,ww] = modelBurst("step", H, Uw, U, 30);
u = ww/(mdl.speed*cos(alfa));
[m, ~] = lsim(Mss, u, t);

figure(2);
subplot(2, 1, 1);
plot(t, u, "LineWidth", 2);
title('Rafaga (Escalon)');
xlabel('t [s]');
ylabel('u [m/s]');
grid on;
hold on;

subplot(2, 1, 2);
plot(t, m + M0, "LineWidth", 2);
xlabel('t [s]');
ylabel('M [Nm]');
title('Respuesta (M_{raiz})');
grid on;
hold on;