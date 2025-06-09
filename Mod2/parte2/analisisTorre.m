close all; clear;

addpath("./functions")

[Mt, Kt, xt] = femTorre();

N = length(xt);

% empotramos la torre
Mt = Mt(3:end, 3:end);
Kt = Kt(3:end, 3:end);

Kt_ = Mt\Kt;
[Vt,Dt] = eig(Kt_, 'vector');
[Dt_s, dtIdx] = sort(Dt);
Wt = sqrt(Dt_s)/(2*pi); %frecuencias de los modos en hz

figure(3)
for h = 1:6
  subplot(3,2,h)
  v = [0;Vt((1:2:2*(N-1)),dtIdx(h))];
  plot(xt, -(1/max(abs(v)))*v, "LineWidth", 2)
  title(["Modo ", num2str(h) " (" num2str(Wt(h)) "Hz)"])

  grid on;
end

