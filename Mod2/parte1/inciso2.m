clear; close all;
pkg load control;

% condición de vuelo
h   = 1000; % m
V   = 140 ; % m/s

addpath('../../Mod1/mdl');
mdl = A320_build_model(1, h, V/1.852*3.6, 0, false, false, 0, 0);
g = mdl.g;
rho = mdl.rho;
U = mdl.speed;
W = mdl.mass;
q = mdl.q;
alfa = mdl.alfa*pi/180;

addpath('../../Mod1/wing');
N = 100;
[M,K,x] = modelwing(N,0);
% cuerda
% area 
% l
% S = sum(areas)

AreaT = 0;

CL = W/2*g/(q*AreaT*alfa);




% TODO: revisar matriz C
C = zeros(1, 2*N);
C(3) = K(3,2);
C(4) = K(4,2);

% TODO: definir matriz F
F = zeros(2*N, 1);

% Empotramos la raiz del ala
M = M(3:end, 3:end);
K = K(3:end, 3:end);
C = C(3:end);
F = F(3:end);

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
P_a = P(idx_a, :); %TODO: entender bien esto
P_b = P(idx_b, :); %TODO: entender bien esto
C_a = C(:, idx_a); %TODO: entender bien esto
C_b = C(:, idx_b); %TODO: entender bien esto

chi = 0.02;

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
[t,u] = modelBurst("1-cos", H, Uw, U, 2.5);
[y, ~] = lsim(Mss, u, t);

figure(2);
subplot(2, 1, 1);
plot(t, u, "LineWidth", 2);
title('Rafaga (1-cos)');
xlabel('t [s]');
ylabel('u [m/s]');
hold on;

subplot(2, 1, 2);
plot(t, y, "LineWidth", 2);
xlabel('t [s]');
ylabel('y [m/s^2]');
title('Respuesta');
hold on;
