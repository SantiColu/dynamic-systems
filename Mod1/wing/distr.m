clear; close all
L = 1;
Rho = 2700;
b = 0.1;
h = 0.1;
A = b*h;
E = 70000*1000;
I = b*h^3/12;

w = 0:2000;
b = nthroot(Rho*A*w.^2/(E*I), 4);

bw = @(W) nthroot(Rho*A*W^2/(E*I), 4);
G = arrayfun(@(W) (cosh(bw(W)*L) + cos(bw(W)*L)).^2 - (sinh(bw(W)*L)-sin(bw(W)*L)).*(sinh(bw(W)*L)+sin(bw(W)*L)), w);
% Gf = @(x) (cosh(x*L) + cos(x*L)).^2 - (sinh(x*L)-sin(x*L)).*(sinh(x*L)+sin(x*L));

bidx = find(G(1:end-1) .* G(2:end) < 0);
bc = b(bidx);
wc = sqrt(((bc.^4)*(E*I))/(Rho*A));
f = wc/(2*pi);

figure(1);
plot(w/(2*pi),G,"LineWidth", 2)
hold on;
for i = 1:length(f)
  fi = f(i);
  plot(fi,0,".", "color", "#2d005b", "MarkerSize", 20)
  if mod(i,2) == 0
    text(fi+5, -10, sprintf("%.2fHz", fi), "FontSize", 14, "fontweight", "bold","VerticalAlignment", "bottom", "HorizontalAlignment", "center")
  else
    text(fi+5, +5, sprintf("%.2fHz", fi), "FontSize", 14, "fontweight", "bold", "VerticalAlignment", "baseline", "HorizontalAlignment", "center")
  end
endfor
grid on;
xlabel('frecuencia [Hz]')
ylabel('G(\beta)')
ylim([-100 100]);

x = linspace(0, L, 100);
length(x)
def = zeros(length(x), 7);

figure(2);
for i = 1:7
  fi = f(i)
  wi = wc(i);
  bet = nthroot(Rho*A*wi^2/(E*I), 4);
  m = -(cosh(bet*L)+cos(bet*L))/(sinh(bet*L)+sin(bet*L));
  B = m;
  y = cosh(bet*x) + m*sinh(bet*x) - cos(bet*x) - m*sin(bet*x);
  y = y/(max(abs(y)));
  if i <= 3
    subplot(2,2,i)
  else
    if i == 4 || i == 5
      subplot(4,4, 7 + i )
    elseif i == 6 || i == 7
      subplot(4,4, 9 + i )
    end
  end
  def(:,i) = y;
  plot(x, y,"LineWidth", 2)
  title(["Modo ", num2str(i) " (" num2str(fi) "Hz)"])
  hold on;
  grid on;
endfor


f = f(1:7);
save("distribuidos.mat", "def", "f");