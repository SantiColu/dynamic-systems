function [t, v] = modelWind(h0, V0, R, w)
  t = 0:0.01:60;
  v = V0^2 * ((30 + 0.75*R*sin(w*t)) / h0).^0.8;
end