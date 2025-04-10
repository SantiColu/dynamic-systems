function [t, x] = run_flight_sim(h, V, gam, slp, d)
  load('./A320_data.mat'); % carga la estructura 'data'
  old_path = addpath('./core');
  switch data.deriv 
      case 'rad'
          data.fderiv = 1;
      case 'deg'
          data.fderiv = 180/pi;
      otherwise
          data.fderiv = 1;
  end
  
  % Como caso particular elegimos como masa el máximo para el aterrizaje
  data.M  = data.mass.max_landing;
  data.dcg = [0 0 0]';
  % El tensor de inercia está estimado para la máxima masa de despegue
  data.J  = data.mass.inertia * (data.M / data.mass.max_takeoff);
  data.Ji = inv(data.J);
  
  %% Condicion de equilibrio
  % densidad según el modelo de atmósfera standard
  rho = densidad(h); 
  % presión dinámica x superficie de referencia
  QS  = 0.5 * rho * V^2 * data.sref;
  % peso
  W   = data.M * 9.81;
  
  % fst(1) = W;
  % fst(2) = gam;
  % fst(3) = QS;
  % % fst(4) = data;
  % fst(5) = 0;
  % % fst(6) = data.flap;
  % % data.flap

  % fst

  trim = level_flight(W, gam, QS, data, 0, data.flap);
  T = trim.T;
  AoA = trim.alfa;
  
  % trimado (empuje y elevador)
  data.T   = [T; 0 ; 0];
  data.cme = -interp1(data.alpha, data.cm, AoA);
  
  % ángulos iniciales [⁰] 
  AoA = AoA + d;    % α [⁰] 
  pch = AoA + gam; % θ [⁰] 
  
  % Vector de estados (6GL)
  % x = { u v w p q r e1 e2 e3 h }
  xo = zeros(10,1);
  
  xo(1)  = V * cosd(AoA);  % u
  xo(2)  = V * sind(slp);  % v
  xo(3)  = V * sind(AoA);  % w
  xo(8)  = pch*pi/180;     % θ [rad]
  xo(10) = h;              % altura [m] 
  
  % usamos ode23, que es menos preciso pero computacionalmente más liviano
  [t, x] = ode23s(@(t,y) xdot(t, y, data) , [0 600], xo); % 
  path(old_path);
end