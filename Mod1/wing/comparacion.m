load("concentrados.mat")
fc = W;
defc = def;

load("distribuidos.mat")
fd = f;
defd = def;

for i = 1:7
  if i <= 3
    subplot(2,2,i)
  else
    if i == 4 || i == 5
      subplot(4,4, 7 + i )
    elseif i == 6 || i == 7
      subplot(4,4, 9 + i )
    end
  end
  yc = defc(:,i);
  yd = defd(:,i);

  frc = fc(i);
  frd = fd(i);
  plot(x, yc, ";concentrados;", "color", "#ff004c", "LineWidth", 2)
  hold on;
  plot(x, yd, ";distribuidos;", "color", "#2443a8", "LineWidth", 2)
  title(["Modo ", num2str(i) " ( f_c = " sprintf("%.2f",frc) " Hz  f_d = " sprintf("%.2f",frd) " Hz   Error = " sprintf("%.1f", abs(frc-frd)/frd*100) "%)"])
  hold on;
  grid on;
endfor