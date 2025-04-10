% x = { u v w p q r e1 e2 e3 h }
% e1: rolido 
% e2: cabeceo 
% e3: rumbo

% 'global' permite declarar variables fuera de la función 
%  que sean accesibles dentro de la función (ver params.m)
function dx = xdot(t, x, data)

    vb = x(1:3); % velocidad lineal
    w  = x(4:6); % velocidad angular
    e  = x(7:9); % euler
    h  = x(10);
    v  = vb + vect_prod(w, data.dcg);
    
    rho = densidad(h);
    
    V2 = v' * v;          % modulo al cuadrado
    Q  = 0.5 * rho * V2;  % presion dinamica
    V  = sqrt(V2);   
    
    a  = atand(v(3)/v(1)); % angulo de ataque
    b  = asind(v(2)/V);    % angulo de deslizamiento
    Ca = stab_to_body(a);
    Ce = dcm(e);
    
    % velocidad angular adimensional 
    p_ = data.blref/(2*V);
    q_ = data.cbar /(2*V);
    w_ = Ca' * w .* [p_ ; q_ ; p_] * data.fderiv; 
    
    fg = Ce * [0 0 9.81]';  % aceleración gravitatoria
    Fa = aero_f(Q, a, b, w_, data); % fuerza aerodinámica
    Fa = Ca * Fa;
    dv = (Fa + data.T)/data.M + fg - vect_prod(w, v);
    
    % da/dt
    tga = v(3)/v(1);
    ad  = 1/(1+tga^2)*(dv(3)/v(1) - dv(1)*v(3)/v(1)^2);
    ad_ = ad * q_;
    
    % momento aerodinámico 
    Ma = aero_m(Q, a, b, ad_, w_, data);
    Ma = Ca * Ma;
    dw = data.Ji*(Ma + vect_prod(Fa, data.dcg) - vect_prod(w, data.J*w));
    
    % derivada de los ángulos de Euler
    de = euler_rate(e, w);
    
    % velocidad vertical
    hdot = -[0 0 1] * Ce'* vb;    
    
    % dx/dt
    dx = [dv ; dw ; de ; hdot]; % .* [1 0 1 0 1 0 0 1 0 1]';
end

% Calcular el producto vectorial entre dos vectores u y v
% lo hacemos usando una matriz anti-simétrica
function w = vect_prod(u, v)
    % calcular el producto vectorial w = u x v
    w = [0    -u(3)  u(2) 
         u(3)  0    -u(1)
        -u(2)  u(1)  0] * v;
end

% calcular fuerzas y momentos aerodinamicos (wind axes)
% Q : presion dinamica 
% a : alfa
% b : beta
% ad: da/dt adimensional
% w : velocidad angular adimensional
function F = aero_f(Q, a, b, w, data)
   
    % calcular lo que corresponda
    cd = interp1(data.alpha, data.cd , a);
    
    cy = interp1(data.alpha, data.cyb, a) * b ...
       + interp1(data.alpha, data.cyp, a) * w(1);
    
    cl = interp1(data.alpha, data.cl , a); ...
    
    F  = Q * data.sref * [-cd ; cy ; -cl];  

end
function M = aero_m(Q, a, b, ad, w, data)

    cl = (interp1(data.alpha, data.clb, a) * b    ...
       +  interp1(data.alpha, data.clp, a) * w(1) ...
       +  interp1(data.alpha, data.clr, a) * w(3)) * data.blref;
   
    cm = (interp1(data.alpha, data.cm , a)        ...
       +  interp1(data.alpha, data.cmad,a) * ad   ...
       +  interp1(data.alpha, data.cmq, a) * w(2) ...
       +  data.cme) * data.cbar;
   
    cn = (interp1(data.alpha, data.cnb, a) * b    ...
       +  interp1(data.alpha, data.cnp, a) * w(1) ...
       +  interp1(data.alpha, data.cnr, a) * w(3)) * data.blref;
   
    M = Q * data.sref * [cl ; cm ; cn];

end

% Calcular derivadas de los angulos de euler
% e: angulos de euler
% w: velocidad de giro
function edot = euler_rate(e, w)
    sn = sin(e(1));
    cs = cos(e(1));
    sc = 1/cos(2);
    tg = sin(2)*sc; 
    C = [1   sn*tg  cs*tg
         0   cs    -sn
         0   sn*sc  cs*sc];
    edot = C * w; % esto es provisorio
end

% Calcular matriz de rotacion
% e: angulos de euler (Tait-Bryan: roll, pitch, yaw)
function C = dcm(e)
    s = sin([e(1), e(2), e(3)]);
    c = cos([e(1), e(2), e(3)]);
    C = [c(2)*c(3)                 c(2)*s(3)                -s(2) 
         s(1)*s(2)*c(3)-c(1)*s(3)  s(1)*s(2)*s(3)+c(1)*c(3)  s(1)*c(2) 
         c(1)*s(2)*c(3)+s(1)*s(3)  c(1)*s(2)*s(3)-s(1)*c(3)  c(1)*c(2)]; 
end

% Transformación wind -> body
function C = stab_to_body(a)
    s = sind(a);
    c = cosd(a);
    C = [ c  0 -s
          0  1  0
          s  0  c];
end

% Modelo de atmósfera estandar
function rho = densidad(H)
    kapha  = -0.000022558; %/m
    rho_sl = 1.225; % a nivel del mar kg/m3
    rho    = rho_sl*((1+kapha*H)^(4.2561)); % a la altura dada kg/m3
end

