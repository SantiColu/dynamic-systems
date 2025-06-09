function [M,K,x] = femTorre()
  f1 = 1.67; % frecuencia del primer modo [Hz]
  b1 = 1.875; % cte (de parametros distribuidos)
  Lt = 24.4; % longitud total de la torre [m]
  mt = 3.2 * 1000; % masa total de la torre [kg]
  rho = mt / Lt; % densidad lineal de masa [kg/m]

  % cantidad de nodos
  N = 20;

  % cantidad de elementos
  n = N-1;

  % coordenadas de los nodos
  x = linspace(0, Lt, N);

  % longitud de cada elemento
  le = Lt/n;

  M = zeros(2*N, 2*N);
  K = zeros(2*N, 2*N);

  ej=(2*pi*f1/b1^2)^2*rho*Lt^4;
  for i = 1:n
    Me = (rho*le)/420 * [156     22*le     54       -13*le
                        22*le    4*le^2    13*le    -3*le^2
                        54       13*le     156      -22*le
                        -13*le   -3*le^2   -22*le   4*le^2];

    Ke = (ej)/le^3 * [12      6*le     -12      6*le
                      6*le    4*le^2   -6*le    2*le^2
                      -12    -6*le     12       -6*le
                      6*le    2*le^2   -6*le    4*le^2];

    idx = 2*(i-1) + (1:4);
    M(idx, idx) = M(idx, idx) + Me;
    K(idx, idx) = K(idx, idx) + Ke;
  endfor
end
