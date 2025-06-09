close all; clear;

addpath("./functions")

[Mp, Kp, rp] = femPala();
N = length(rp);

% empotramos la pala
Mp = Mp(3:end, 3:end);
Kp = Kp(3:end, 3:end);

Kp_ = Mp\Kp;
[Vp,Dp] = eig(Kp_, 'vector');
[Dp_s, dpIdx] = sort(Dp);
Wp = sqrt(Dp_s)/(2*pi); %frecuencias de los modos en hz

figure(3)
for h = 1:6
  subplot(3,2,h)
  v = [0;Vp((1:2:2*(N-1)),dpIdx(h))];
  plot(rp, -(1/max(abs(v)))*v, "LineWidth", 2)
  title(["Modo ", num2str(h) " (" num2str(Wp(h)) "Hz)"])

  grid on;
end

