function [t, u] = modelBurst(typ, H, Uw, U, periods = 2.5)
  if strcmp(typ, '1-cos')
    t = 0:0.005:periods*(2*H/U);
    tm = zeros(size(t));
    tm(t < 2*H/U) = t(t < 2*H/U);
    u = Uw * sin(pi*U*tm/(2*H)).^2;
    return
  end


  % TODO: quieren escalon o salto cte?
  if strcmp(typ, 'step')
    t = 0:0.005:periods*(2*H/U);
    u = Uw * ones(size(t)) .* (t > 2*H/U & t < 2*2*H/U);
    return
  end
end