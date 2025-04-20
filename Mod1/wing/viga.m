clear;

N = 5;
L = 1;

b = 0.1;
h = 0.1;
A = b*h;
Rho = 2700;

m = A*Rho*L;

E = 70000*1000;
J = b*h^3/12;

me = m/(N-1);
Le = L/(N-1);
Ie = me*Le^2/12;

M = zeros(2*N);
for i = 1:N
  if i == 1 || i == N
    M(2*i - 1, 2*i - 1) = me/2;
    M(2*i, 2*i) = me/2*(Le/2+h)^2/3; %modificar
  else
    M(2*i - 1, 2*i - 1) = me;
    M(2*i, 2*i) = me*(Le+h)^2/12;
  end
end

Ke = [12/(Le^3)    6/(Le^2)    -12/(Le^3)   6/(Le^2);
      6/(Le^2)     4/Le        -6/(Le^2)    2/Le;
      -12/(Le^3)   -6/(Le^2)   12/(Le^3)    -6/(Le^2);
      6/(Le^2)     2/Le        -6/(Le^2)    4/Le];

K = zeros(2*N);
for i = 1:2:2*(N-1)
  rng = [0 1 2 3] + i;
  K(rng, rng) = K(rng, rng) + Ke;
end

M = M(3:end, 3:end);
K = E*J*K(3:end, 3:end);

K_ = M\K;
[V,D] = eig(K_, 'vector');

[D_s, dIdx] = sort(D);
W = sqrt(D_s)/(2*pi); %frecuencias de los modos en hz

x = (0:N-1)' * Le;
length(x)

def = zeros(length(x), 7);

figure(3)
for h = 1:7
  if h <= 3
    subplot(2,2,h)
  else
    if h == 4 || h == 5
      subplot(4,4, 7 + h )
    elseif h == 6 || h == 7
      subplot(4,4, 9 + h )
    end
  end
  v = [0;V((1:2:2*(N-1)),dIdx(h))];
  t = [0;V((2:2:2*(N-1)+1),dIdx(h))];
  if h == 3 || h == 5 || h == 2 || h == 7
    v = -v;
  end
  def(:,h) = (1/max(abs(v)))*v;
  plot(x, def(:,h), "LineWidth", 2)
  title(["Modo ", num2str(h) " (" num2str(W(h)) "Hz)"])

  % plot(x, (1/max(abs(t)))*t, ["r;giros modo " num2str(h) ";"])
  grid on;
end

W = W(1:7);
save("concentrados.mat", "def", "W");