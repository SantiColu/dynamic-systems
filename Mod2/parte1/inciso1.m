clear; clc; close all;
pkg load control;
addpath("./functions")

% condición de vuelo
h   = 1000; % m
V   = 140 ; % m/s

addpath('../../Mod1/mdl');
mdl = A320_build_model(1, h, V/1.852*3.6, 0, false, false, 0, 0);
U = mdl.speed * cosd(mdl.alfa);

A = mdl.lng.A;
B = mdl.lng.B;
C = A(2,:);
C(3) = 0;
D = B(2,:);
M = ss(A, B, C, D, 'StateName', mdl.lng.StateName, 'InputName', mdl.lng.InputName, 'OutputName', {"FC"});
M_ = M(1,4);
G = tf(M_)


% Para la rafaga (1-cos)
Hr = (30:50:350)*0.3048; % m
Uw = 10; % m/s^2
Ymax = [0,0];
for i = 1:length(Hr)
  H = Hr(i);
  [t,u] = modelBurst("1-cos",H, Uw, U, 5);
  [y, ~] = lsim(M_, u, t);

  figure(1);
  subplot(2, 1, 1);
  plot(t, u, "LineWidth", 2, sprintf(';H = %.1f m;', H));
  title('Rafaga (1-cos)');
  xlabel('t [s]');
  ylabel('u [m/s]');
  grid on;
  hold on;

  subplot(2, 1, 2);
  plot(t, y-mdl.g, "LineWidth", 2,sprintf(';H = %.1f m;', H));
  xlabel('t [s]');
  ylabel('FC [m/s^2]');
  title('Respuesta');
  grid on;
  hold on;

  [ym, im] = max(abs(y-mdl.g));
  if ym > abs(Ymax(2))
    Ymax = [i, y(im)-mdl.g];
  end
endfor

subplot(2, 1, 2);
plot(xlim(), [Ymax(2), Ymax(2)], 'r--;FC_{max};', "LineWidth", 1.5);
text(xlim()(2), Ymax(2)*0.9, sprintf('Ymax = %.2f m/s^2', Ymax(2)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', "color", "red");


% Para el escalon
H = Hr(Ymax(1));
[t,u] = modelBurst("step",H, Uw, U, 15);
[y, ~] = lsim(M_, u, t);

figure(2);
subplot(2, 1, 1);
plot(t, u, "LineWidth", 2);
title('Rafaga (escalon)');
xlabel('t [s]');
ylabel('u [m/s]');
grid on;
hold on;

subplot(2, 1, 2);
plot(t, y-mdl.g, "LineWidth", 2);
xlabel('t [s]');
ylabel('FC [m/s^2]');
title('Respuesta');
grid on;
hold on;

% TODO: revisar metodo 2
