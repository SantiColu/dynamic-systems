clear; close all;

N = 100;
[M,K,x] = modelwing(N,1);

M = M(3:end, 3:end);
K = K(3:end, 3:end);

K_ = M\K;

[V,D] = eig(K_, 'vector');
[D_s, dIdx] = sort(D);
W = sqrt(D_s)/(2*pi); %frecuencias de los modos en hz

figure(3)
for h = 1:6
  subplot(3,2,h)
  v = [0;V((1:2:2*(N-1)),dIdx(h))];
  % t = [0;V((2:2:2*(N-1)+1),dIdx(h))];
  plot(x, -(1/max(abs(v)))*v, "LineWidth", 2)
  title(["Modo ", num2str(h) " (" num2str(W(h)) "Hz)"])

  % plot(x, (1/max(abs(t)))*t, ["r;giros modo " num2str(h) ";"])
  grid on;
end

