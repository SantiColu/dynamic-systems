clear; clc; clf

load('../A320_data.mat'); % carga la estructura 'data'
old_path = addpath('../core');

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
    
    % condición de vuelo
    h   = 5000; % m
    V   = 130 ; % m/s
    gam = 0;    % γ
    slp = 0;    % β
    d   = 1;    % Δα / Δθ [⁰] - apartamiento inicial drespecto del equilibrio
    
    
    %% Condicion de equilibrio
    % densidad según el modelo de atmósfera standard
    rho = densidad(h); 
    % presión dinámica x superficie de referencia
    QS  = 0.5 * rho * V^2 * data.sref;
    % peso
    W   = data.M * 9.81;
    
    trim = level_flight(W, gam, QS, data, 0, data.flap);
    T = trim.T;
    AoA = trim.alfa;
    
    % trimado (empuje y elevador)
    data.T   = [T; 0 ; 0];
    data.cme = -interp1(data.alpha, data.cm, AoA);
    
    % ángulos iniciales [⁰] 
    AoA = AoA + d    % α [⁰] 
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
    
    %%
    clf
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

path(old_path);



