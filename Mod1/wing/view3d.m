% Call model
n = 5;
L = 15.075;
Le = L/(n-1);
x = Le*(0:n-1);

show_rect = true;

airfoil = load("airfoil.dat");
x_airf = airfoil(:, 1);
y_airf = airfoil(:, 2);
x_airf = x_airf / max(x_airf);

% Wing
L = 15.075;
brng = [6, 1.5];
trng = [0.008, 0.005];

figure(1)
clf
hold on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('Modelo 3D del ala')
view(55, 15)
lighting gouraud
camlight headlight

for i = 1:n-1
  x1 = x(i);
  x2 = x(i+1);

  xpos = x1;
  b = interp1([0, L], brng, xpos);
  t = interp1([0, L], trng, xpos);

  xa = x_airf * b;
  ya = y_airf * b;

  h = max(ya) - min(ya);

  % === airfoil ===
  for j = 1:length(xa)-1
    v1 = [x1, xa(j),   ya(j)];
    v2 = [x2, xa(j),   ya(j)];
    v3 = [x2, xa(j+1), ya(j+1)];
    v4 = [x1, xa(j+1), ya(j+1)];
    patch('Vertices', [v1; v2; v3; v4], ...
          'Faces', [1 2 3 4], ...
          'FaceColor', [0.8 0.8 0.8], ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.95);
  end

  % === rectangular approximation===
  if show_rect
    y_rect = [0, b];
    z_rect = [-h/2, h/2];
    for j = 1:4
      switch j
        case 1  % bottom
          y = y_rect;
          z = z_rect(1) * ones(1,2);
        case 2  % top
          y = y_rect;
          z = z_rect(2) * ones(1,2);
        case 3  % left
          y = y_rect(1) * ones(1,2);
          z = z_rect;
        case 4  % right
          y = y_rect(2) * ones(1,2);
          z = z_rect;
      endswitch

      v1 = [x1, y(1), z(1)];
      v2 = [x2, y(1), z(1)];
      v3 = [x2, y(2), z(2)];
      v4 = [x1, y(2), z(2)];

      patch('Vertices', [v1; v2; v3; v4], ...
          'Faces', [1 2 3 4], ...
          'FaceColor', [0.6 0.8 1], ...
          'EdgeColor', [0.2 0.2 0.6], ...
          'LineStyle', '--', ...
          'LineWidth', 0.5);
    end
  end
end


grid on

% Motor 
Mmotor = 2500;
XMotor = 3.5;  % m
L = 15.075;
brng = [6, 1.5];
trng = [0.008, 0.005];

b_motor = interp1([0, L], brng, XMotor);
t_motor = interp1([0, L], trng, XMotor);

airfoil = load("airfoil.dat");
x_airf = airfoil(:, 1);
y_airf = airfoil(:, 2);
x_airf = x_airf - min(x_airf);
x_airf = x_airf / max(x_airf);
xa = x_airf * b_motor;
ya = y_airf * b_motor;
p = 0;
for k = 1:length(xa)-1
  p += sqrt((xa(k+1) - xa(k))^2 + (ya(k+1) - ya(k))^2);
end
h_motor = p/2 - b_motor;

x_mot = XMotor;
y_mot = 0;
z_mot = 0;

[xx, yy, zz] = sphere(10);
r = min(b_motor, h_motor) * 0.7;
surf(xx*r + x_mot, yy*r + y_mot, zz*r + z_mot, ...
     'FaceColor', [1 0 0], ...
     'EdgeColor', 'none');