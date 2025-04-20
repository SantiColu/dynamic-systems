function [M,K,x] = modelwing(n, fig)
  L = 15.075; %m
  RhoAl = 2700; %kg/m^3
  RhoJA1 = 800; %kg/m^3
  EJalfa = 750000000;
  Mmotor = 2500; %kg
  XMotor = 3.5; %m

  Le = L/(n-1);

  brng = [6, 1.5]; %m
  trng = [0.008, 0.005]; %m (es un poco mas para añadir largueros y costillas (?))
  x = Le*(0:n-1);

  ixmot = round(XMotor/Le) + 1;

  airfoil = load("airfoil.dat");
  x_airf = airfoil(:, 1);
  y_airf = airfoil(:, 2);

  if fig != 0
    figure(fig)
    plot(x_airf, y_airf, "LineWidth", 2)
    hold on
    axis equal
    grid on
  end

  Ke = [12/(Le^3)    6/(Le^2)    -12/(Le^3)   6/(Le^2);
        6/(Le^2)     4/Le        -6/(Le^2)    2/Le;
        -12/(Le^3)   -6/(Le^2)   12/(Le^3)    -6/(Le^2);
        6/(Le^2)     2/Le        -6/(Le^2)    4/Le];
  J = zeros(1,n);

  M = zeros(2*n);
  K = zeros(2*n);

  masa = 0;
  m = zeros(n-1, 1);
  I = zeros(n-1, 1);
  for i = 1:n
    Lei = Le;
    if i == 1 || i == n
      Lei = Le/2;
    end
    xpos = (i-1)*Le;
    b = interp1([0, L], brng, xpos);
    t = interp1([0, L], trng, xpos);
    nairfx = x_airf*b;
    nairfy = y_airf*b;

    %calculate the airfoil perimeter
    p=0;
    for k = 1:length(nairfx)-1
      p = p + sqrt((nairfx(k+1) - nairfx(k))^2 + (nairfy(k+1) - nairfy(k))^2);
    end

    h = max(nairfy) - min(nairfy);
    m(i) = p*t*Lei*RhoAl; %add fuel mass
    if i == ixmot
      m(i) = m(i) + Mmotor;
    end
    I(i) = 1/12*m(i)*Lei^2;

    M(2*i - 1, 2*i - 1) =  m(i);
    M(2*i, 2*i) = I(i);
    masa = masa + m(i);

    J(i) = (b*h^3)/12;
  endfor

  c = 1;
  for i = 1:2:2*(n-1)
    Jp = (J(c)+J(c+1))/2;
    rng = [0 1 2 3] + i;
    K(rng, rng) = K(rng, rng) + EJalfa*Jp*Ke;
    c = c + 1;
  end

  figure(fig+1)
  subplot(311)
  plot(x, m,";Masa;", "LineWidth", 2)
  hold on;
  plot([XMotor, XMotor], ylim(), "r--;motor;", "LineWidth", 1.5)
  xlabel('x [m]')
  ylabel('Masa del elemento [kg]')
  % plot([x(ixmot), x(ixmot)], ylim(), "m--;motor x;", "LineWidth", 1.5)
  subplot(312)
  plot(x, I,";I;", "LineWidth", 1.5)
  hold on;
  plot([XMotor, XMotor], ylim(), "r--;motor;", "LineWidth", 1.5)
  subplot(313)
  plot(x, J,";J;", "LineWidth", 2)
  hold on;
  plot([XMotor, XMotor], ylim(), "r--;motor;", "LineWidth", 1.5)
  xlabel('x [m]')
  ylabel('J [m^4]')
end
