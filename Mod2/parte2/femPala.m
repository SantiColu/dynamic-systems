function [M,K,r,c] = femPala()
  load("AOC_15_50.mat")

  % longitud de cada elemento
  L = diff(r);

  % cantidad de elementos
  n = length(L);

  % cantidad de nodos
  N = length(r);

  M = zeros(2*N, 2*N);
  K = zeros(2*N, 2*N);

  for i = 1:n
    le = L(i);
    ej = EJ(i) * 10^6;
    rho = roA(i);

    Me = (rho*le)/420 * [156     22*le     54       -13*le
                        22*le    4*le^2    13*le    -3*le^2
                        54       13*le     156      -22*le
                        -13*le   -3*le^2   -22*le   4*le^2];

    Ke = (ej)/le^3 * [12     6*le     -12      6*le
                      6*le    4*le^2   -6*le    2*le^2
                      -12    -6*le     12       -6*le
                      6*le    2*le^2   -6*le    4*le^2];

    idx = 2*(i-1) + (1:4);
    M(idx, idx) = M(idx, idx) + Me;
    K(idx, idx) = K(idx, idx) + Ke;
  endfor
end
