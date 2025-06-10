function [v_aprox, m, a0] = fftWind(t, v, pl)
  if nargin < 3
    pl = false;  % Por defecto, graficar los coeficientes
  end

  N = length(t);
  dt = t(2) - t(1);     % Paso de tiempo
  fs = 1/dt;            % Frecuencia de muestreo
  f = fs * (1:(N/2 - 1)) / N;  % Frecuencias en Hz (hasta Nyquist)

  % FFT y coeficientes
  U = fft(v)/N;
  a0 = real(U(1));
  an = 2*real(U(2:N/2));
  bn = -2*imag(U(2:N/2));

  if pl
    % Graficar an y bn con eje de frecuencias
    figure;
    subplot(2,1,1);
    bar(f, an(1:length(f)), 'b');
    xlabel('Frecuencia [Hz]'); ylabel('a_n');
    title('Coeficientes coseno a_n'); grid on;
    xlim([0 4]);

    subplot(2,1,2);
    bar(f, bn(1:length(f)), 'r');
    xlabel('Frecuencia [Hz]'); ylabel('b_n');
    title('Coeficientes seno b_n'); grid on;
    xlim([0 4]);
  end

  % Aproximación con el mayor bn
  [m, i] = max(bn);
  f_i = f(i);  % frecuencia dominante
  v_aprox = a0 + m * sin(2*pi*f_i*t);
end