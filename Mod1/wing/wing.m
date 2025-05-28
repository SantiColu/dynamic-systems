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


% -------- Modelo modo 1 --------

modo1 = V(:, dIdx(1));
modo1 = modo1 / max(abs(modo1));

m_eq = (modo1' * M * modo1) ;
k_eq = (modo1' * K * modo1) ;
f_eq = sqrt(k_eq / m_eq) / (2*pi);

fprintf('\n-------- Modelo modo 1 --------\n');
fprintf('  m_eq:  %.2f kg\n', m_eq);
fprintf('  k_eq:  %.2e N/m\n', k_eq);
fprintf('  f_eq:  %.2f Hz\n', f_eq);