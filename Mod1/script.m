clear; clc; close all;
addpath('./sim');

force_sim_run = false;

% condición de vuelo
h   = 5000; % m
V   = 130 ; % m/s
gam = 0;    % γ
slp = 1;    % β
d   = 0;    % Δα / Δθ [⁰]
phi = 1; % [⁰] roll

simName = sprintf("output/flight_sim_h%.0f_v%.1f_gam%.1f_slp%.1f_d%.1f_ps%.1f_ph%.1f.mat", h, V, gam, slp, d, 0, phi);

% ----------------------------  MODELO REAL ---------------------------- 
  if exist(simName, 'file')
    if force_sim_run
      disp('simulation found, but forced to run again...');
      load_real = false;
    else
      disp('simulation found, loading...');
      load_real = true;
    end
  else
    disp('simulation not found, running sim...');
    load_real = false;
  end

  if load_real
    load(simName)
  else
    [t, x] = run_flight_sim(h, V, gam, slp, d, psi*pi/180, phi*pi/180);
    save(simName, 'x', 't');
  end
  figure(1);
  subplot(421); plot(t, x(:,1)); grid on; ylabel('u [m/s]');
  subplot(423); plot(t, atand(x(:,3)./x(:,1))); grid on; ylabel('\alpha [⁰]'); %ylim([0 10])
  subplot(425); plot(t, x(:,8)*180/pi); grid on; ylabel('\theta [⁰]');
  subplot(427); plot(t, x(:,10)); grid on; ylabel('h [m]');

  xlabel('tiempo');

  subplot(422); plot(t, x(:,[4 6])*180/pi); grid on; ylabel('p/r [⁰/s]');
  subplot(424); plot(t, atand(x(:,2)./x(:,1))); grid on; ylabel('\beta [⁰]');
  subplot(426); plot(t, x(:,7)*180/pi); grid on; ylabel('\phi [⁰]');
  subplot(428); plot(t, x(:,9)*180/pi); grid on; ylabel('\psi [⁰]');

  xlabel('tiempo');
% ----------------------------  MODELO REAL ---------------------------- 



% ---------------------------- MODELO LINEAL ---------------------------
  addpath('./mdl');
  mdl = A320_build_model(1, h, V/1.852*3.6, 0, false, false, gam, slp);

  % lng:
  X0_lng = [0; d*130*pi/180; 0; (d)*pi/180];
  A_lng = mdl.lng.A;
  [V_lng,D_lng] = eig(A_lng);
  A_lng_ = V_lng\(A_lng*V_lng);
  Z0_lng = V_lng\X0_lng;
  Z_lng = zeros(4,length(t));
  X_lng = zeros(4,length(t));

  % lng:
  X0_lat = [slp*pi/180; 0; 0; phi*pi/180];
  A_lat = mdl.lat.A;
  [V_lat,D_lat] = eig(A_lat);
  Z0_lat = V_lat\X0_lat;
  Z_lat = zeros(4,length(t));
  X_lat = zeros(4,length(t));

  for i = 1:4
    Z_lat(i, :) = Z0_lat(i) * e.^(D_lat(i,i)*t);
    Z_lng(i, :) = Z0_lng(i) * e.^(D_lng(i,i)*t);
  endfor

  X_lat = V_lat*Z_lat;
  X_lng = V_lng*Z_lng;

  mdl_res.u = X_lng(1, :) + mdl.speed * cosd(mdl.alfa);
  mdl_res.w = X_lng(2, :) + mdl.speed * sind(mdl.alfa);
  mdl_res.alfa = atand(mdl_res.w./mdl_res.u);
  mdl_res.pitch = X_lng(4, :) + mdl.pitch*pi/180;

  mdl_res.h(1) = h;
  for i = 2:length(t)
    mdl_res.h(i) = mdl_res.h(i-1) + X_lng(2, i)*(t(i)-t(i-1)); % TODO: REVISAR
  endfor

  figure(2);
  subplot(421); plot(t, mdl_res.u); grid on; title("u")
  subplot(423); plot(t, mdl_res.alfa); grid on; title("alfa")
  subplot(425); plot(t, mdl_res.pitch*180/pi); grid on; title("pitch (tita)")
  subplot(427); plot(t, mdl_res.h); grid on; title("h")

  subplot(422); plot(t, X_lat(1, :)); grid on; title("beta")
  subplot(424); plot(t, X_lat(2, :)); grid on; title("p")
  subplot(426); plot(t, X_lat(3, :)); grid on; title("r")
  subplot(428); plot(t, X_lat(4, :)); grid on; title("roll")
% ---------------------------- MODELO LINEAL ---------------------------


% ---------------------------- COMPARACION ---------------------------
  figure(3);
  subplot(421);
    plot(t, x(:,1), "color", "#2443a8", ";sim;", "LineWidth", 2);
    hold on;
    plot(t, mdl_res.u, "color", "#ff004c", ";lineal;", "LineWidth", 2);
    grid on; 
    ylabel('u [m/s]');

  subplot(423); 
    plot(t, atand(x(:,3)./x(:,1)), "color", "#2443a8", ";sim;", "LineWidth", 2);
    hold on;
    plot(t, mdl_res.alfa, "color", "#ff004c", ";lineal;", "LineWidth", 2);
    grid on; 
    ylabel('\alpha [⁰]'); 

  subplot(425); 
    plot(t, x(:,8)*180/pi, "color", "#2443a8", ";sim;", "LineWidth", 2);
    hold on;
    plot(t, mdl_res.pitch*180/pi, "color", "#ff004c", ";lineal;", "LineWidth", 2);
    grid on; 
    ylabel('\theta [⁰]');

  subplot(427);
    plot(t, x(:,3), "color", "#2443a8", ";sim;", "LineWidth", 2);
    hold on;
    plot(t, mdl_res.w, "color", "#ff004c", ";lineal;", "LineWidth", 2);
    grid on;
    ylabel('w [m/s]');



  subplot(422);
    plot(t, x(:,4)*180/pi, "color", "#2443a8", ";sim;", "LineWidth", 1.5);
    hold on;
    plot(t, X_lat(2, :)*180/pi, "color", "#ff004c", ";lineal;", "LineWidth", 1.5);
    grid on; 
    ylabel('p [⁰/s]');

  subplot(424);
    hold on;
    plot(t, x(:,6)*180/pi, "color", "#2443a8", ";sim;", "LineWidth", 1.5);
    plot(t, X_lat(3, :)*180/pi, "color", "#ff004c", ";lineal;", "LineWidth", 1.5);
    grid on; 
    ylabel('r [⁰/s]');

  subplot(426); 
    plot(t, atand(x(:,2)./x(:,1)), "color", "#2443a8", ";sim;", "LineWidth", 1.5);
    hold on;
    plot(t, X_lat(1, :)*180/pi, "color", "#ff004c", ";lineal;", "LineWidth", 1.5);
    grid on; 
    ylabel('\beta [⁰]'); 

  subplot(428); 
    plot(t, x(:,7)*180/pi, "color", "#2443a8", ";sim;", "LineWidth", 1.5);
    hold on;
    plot(t, X_lat(4, :)*180/pi, "color", "#ff004c", ";lineal;", "LineWidth", 1.5);
    grid on; 
    ylabel('\phi [⁰]');

% ---------------------------- ANALISIS ---------------------------
original_t = t;
function [max_indices, min_indices] = encontrar_extremos(y, include_first)
  if nargin < 2
    include_first = false;
  end

  max_indices = [];
  min_indices = [];

  if include_first
    if y(1) > y(2) 
      max_indices(1) = 1;
    elseif y(1) < y(2)
      min_indices(1) = 1;
    end
  end

  for i = 2:length(y)-1
    if y(i) > y(i-1) && y(i) > y(i+1)
      max_indices(end+1) = i;
    elseif y(i) < y(i-1) && y(i) < y(i+1)
      min_indices(end+1) = i;
    end
  end
end

function plot_extremos(fig, t, y, max_indices, min_indices, leg, xa, ya, lineal)
  figure(fig);
  h(1) = plot(t,y, "color", "#ff004c", [";" leg ";"], "LineWidth", 1.5);
  hold on;
  if nargin >= 9
    h(2) = plot(t,lineal, "color", "#2443a8", [";lineal;"], "LineWidth", 1.5);
  end

  grid on;
  ylims = ylim();

  for i = max_indices
      plot([t(i), t(i)], ylims, 'b--'); 
  endfor

  for i = min_indices
      plot([t(i), t(i)], ylims, 'm--'); 
  endfor
  legend(h)
  xlabel(xa);
  ylabel(ya);
end


%  --- MODO LENTO --- 

% frecuencia
omega = imag(D_lng(3,3));
f_lin = omega/(2*pi);

% amortiguamiento
amort_lin = real(D_lng(3,3));

% utilizamos u (sim para conseguir la frecuencia y el factor de amortiguamiento)
y = x(:,1)(100:end);
t = t(100:end);

[max_indices, min_indices] = encontrar_extremos(y);
plot_extremos(4, t, y, max_indices, min_indices, "u", "tiempo [s]", "u [m/s]");

% obtenemos los valores y tiempos de los maxs/mins
% nos sacamos los primeros 2 porque tienen tambien el modo lento
y_max = y(max_indices)(3:end);
t_max = t(max_indices)(3:end);
y_min = y(min_indices)(3:end);
t_min = t(min_indices)(3:end);

% sacamos el promedio del tiempo entre maximos y minimos
Tmi = t_min(5)-t_min(4);
Tma = t_max(6)-t_max(5);
T = ((Tma+Tmi)/2);
f = 1/T;

% para el amortiguamiento, bajamos la curva para que converja a 0
m = (y_max(1)+y_min(1))/2;
y_max = y_max - m;

% hacemos un fit exponencial y obtenemos las constantes
lny = log(y_max);
p = polyfit(t_max, lny, 1);
amort = p(1);
a = exp(p(2));
fitted_y = a * exp(amort * t);

figure(5);
plot(t,y-m, ";sim;", "LineWidth", 1.5);
hold on;
plot(t,mdl_res.u(100:end)-m, ";lin;", "LineWidth", 1.5);
plot(t,fitted_y, ";fitted_y;", "LineWidth", 1.5);
xlabel("tiempo [s]");
ylabel("u [m/s]");
grid on;

% calculamos los errores
error_f = (abs(f - f_lin)/f)*100;
error_a = (abs(amort_lin - amort)/abs(amort))*100;

disp("\n--- MODO LENTO ---")
disp(sprintf("f_real: %f,\t f_lin: %f,\t error: %.1f%%", f, f_lin, error_f));
disp(sprintf("amort_real: %f,\t amort_lin: %f,\t error: %.1f%%", amort, amort_lin, error_a));
disp("------------------\n")

t = original_t;
%  --- MODO RAPIDO --- 

% frecuencia
omega = imag(D_lng(1,1));
f_lin = omega/(2*pi);

% amortiguamiento
amort_lin = real(D_lng(1,1));


% utilizamos alfa (sim para conseguir la frecuencia y el factor de amortiguamiento)
y = atand(x(:,3)./x(:,1))(1:70);
t = t(1:70);

[max_indices, min_indices] = encontrar_extremos(y, 1);
plot_extremos(6, t, y, max_indices, min_indices, "alfa", "tiempo [s]", '\alpha [⁰]');

% obtenemos los valores y tiempos de los maxs/mins
y_max = y(max_indices);
t_max = t(max_indices);
y_min = y(min_indices);
t_min = t(min_indices);

% sacamos el promedio del tiempo entre maximos y minimos
Tmi = t_min(2)-t_min(1);
Tma = t_max(2)-t_max(1);
T = ((Tma+Tmi)/2);
f = 1/T;

% para el amortiguamiento, bajamos la curva para que converja a 0
m = y(end);
y_max = y_max - m;

% hacemos un fit exponencial y obtenemos las constantes
lny = log(y_max);
p = polyfit(t_max, lny, 1);
amort = p(1);
a = exp(p(2));
fitted_y = a * exp(amort * t);

figure(7);
plot(t,y-m, ";sim;","LineWidth", 1.5);
hold on;
plot(t,mdl_res.alfa(1:70)-m, ";lin;", "LineWidth", 1.5);
plot(t,fitted_y, ";fitted_y;","LineWidth", 1.5);
grid on;
xlabel("tiempo [s]");
ylabel('\alpha [⁰]');

error_f = (abs(f - f_lin)/f)*100;
error_a = (abs(amort_lin - amort)/abs(amort))*100;

disp("\n--- MODO RAPIDO ---")
disp(sprintf("f_real: %f,\t f_lin: %f,\t error: %.1f%%", f, f_lin, error_f));
disp(sprintf("amort_real: %f,\t amort_lin: %f,\t error: %.1f%%", amort, amort_lin, error_a));
disp("---------------------\n")

t = original_t;
%  --- LATERAL --- 

% frecuencia
omega = imag(D_lat(1,1));
f_lin = omega/(2*pi);

% utilizamos beta (sim para conseguir la frecuencia)
y = atand(x(:,2)./x(:,1))(1:300);
t = t(1:300);

[max_indices, min_indices] = encontrar_extremos(y, 1);
plot_extremos(8, t, y, max_indices, min_indices, "beta", "tiempo [s]", '\beta [⁰]', (X_lat(1, :)*180/pi)(1:300));

% obtenemos los valores y tiempos de los maxs/mins
y_max = y(max_indices);
t_max = t(max_indices);
y_min = y(min_indices);
t_min = t(min_indices);

% sacamos el promedio del tiempo entre maximos y minimos
Tmi = t_min(2)-t_min(1);
Tma = t_max(2)-t_max(1);
T = ((Tma+Tmi)/2);
f = 1/T;

error_f = (abs(f - f_lin)/f)*100;
error_a = (abs(amort_lin - amort)/abs(amort))*100;

disp("\n--- LATERAL ---")
disp(sprintf("f_real: %f,\t f_lin: %f,\t error: %.1f%%", f, f_lin, error_f));
disp("---------------------\n")
