clear; close all;

N = 7;
m = 1;
L = 1;
T = 1;

Le = L/(N-1);
Ke = [1 -1; -1 1];

M = eye(N) * m/(N-2);
K = zeros(N);
for i = 1:N-1
  rng = [0 1] + i;
  K(rng, rng) = K(rng, rng) + Ke;
end

K = K(2:end-1, 2:end-1);
M = M(2:end-1, 2:end-1);

K_ = M\K;
[V,D] = eig(K_, 'vector');

x = (0:N-1)' * Le;

for h = 1:3
  v = [0;V(:,h);0];
  plot(x, -v, [";modo " num2str(h) ";"])
  hold on
end

