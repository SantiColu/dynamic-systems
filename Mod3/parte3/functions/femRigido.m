function [G, mag] = femRigido(V0, w)
  load("functions/ASW20_flex_1.mat");
  load("functions/ASW20_rigid_long.mat"); %(Ai, Aa, M)

  la=15/2;
  c=0.7;
  CLa=2*pi;
  dca=2.42-2.29;
  rho=1.204;
  g=9.81;

  gam = 0.5*rho*V0*c*CLa;

  A = Ai+Aa;
  B = M\[0; -gam*la; -gam*dca*la; 0];
  C = A(2,:);
  C(3)=0;
  C(4)=0;
  D = B(2,:);

  M = ss(A, B, C, D);
  G = tf(M);
  mag=bode(G,w);
  mag=squeeze(mag);
end

