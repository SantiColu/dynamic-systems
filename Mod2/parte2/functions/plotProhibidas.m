function plotProhibidas(fig, f, xl, yl, ti)
  % Calcular márgenes ±10%
  bandas = [f(:)*0.90, f(:)*1.10];
  % Ordenar por frecuencia de inicio
  bandas = sortrows(bandas);

  % Unir bandas superpuestas
  bandas_unidas = [];
  actual = bandas(1,:);
  for i = 2:size(bandas, 1)
      if bandas(i,1) <= actual(2)  % Se solapan
          actual(2) = max(actual(2), bandas(i,2));  % Unir
      else
          bandas_unidas = [bandas_unidas; actual];
          actual = bandas(i,:);
      end
  end
  bandas_unidas = [bandas_unidas; actual];  % Añadir última banda

  % Rango total del gráfico
  f_min = 0;
  f_max = max(f)*1.5;

  % Graficar
  hf = figure(fig);
  hold on;
  ylim([0, 1]);
  xlim([f_min, f_max]);

  % Dibujar zonas rojas unidas
  for i = 1:size(bandas_unidas, 1)
      f_low = bandas_unidas(i, 1);
      f_high = bandas_unidas(i, 2);
      fill([f_low f_high f_high f_low], [0 0 1 1], [1 0 0], ...
          'FaceAlpha', 0.3, 'EdgeColor', 'none');
  end

  % Dibujar líneas verticales en frecuencias originales
  for i = 1:length(f)
      plot([f(i) f(i)], [0 1], 'r--', 'LineWidth', 1.2);
  end

  % Etiquetas
  xlabel(xl);
  ylabel(yl);
  title(ti);
  yticks([]);
  grid on;

  print (hf, ['/home/santi/UNLP/SistDinamicos/Mod2/parte2/presentacion/frecuencias/plotProhibidas_', num2str(fig), '.jpg'], "-djpeg");
end