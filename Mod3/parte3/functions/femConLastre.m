function [G, mag] = femConLastre(V0, w)
  load("functions/ASW20_flex_1.mat");
  load("functions/ASW20_rigid_long.mat"); %(Ai, Aa, M)

  I1=trapz(X, Y1);
  I2=trapz(X, Y1.^2);


  la=15/2;
  c=0.7;
  CLa=2*pi;
  dca=2.42-2.29;
  rho=1.204;
  g=9.81;

  %% Modelo semi-ala con lastre
  fa=1.77;                         %Frecuencia natural [Hz]
  wa=2*pi*fa;                      %Frecuencia angular natural [rad/s]
  e=0.01;                          %Factor de amortiguamiento estructural

  ma=65 + 60;                      %Masa semi-ala con lastre [kg]
  rhoA=ma/la;                      %Densidad lineal semi-ala [kg/m]
  meq=rhoA*I2;                     %Masa equivalente semi-ala 
  keq=wa^2*meq;                    %Rigidez equivalente semi-ala 
  deq=2*e*wa*meq;                  %Disipación equivalente semi-ala 


  V0=120*1000/3600;
  gam = 0.5*rho*V0*c*CLa*2;

  A11 = Ai + Aa;
  B1 = M\[0; -gam*la; -gam*dca*la; 0];

  Qze = -gam*I1;
  Qte = gam*dca*I1;
  Qee = -gam*I2;
  A12 = M\[0 0; 0 Qze; 0 Qte; 0 0];

  Qer_w = -gam*I1;
  Qer_q = -gam*dca*I1;
  A21 = [ 0  0          0          0; 
          0  Qer_w/meq  Qer_q/meq  0];

  A22 = [ 0      1;
          -wa^2  -(deq - Qee)/meq];

  Qew = -gam*I1;
  B2 = [0; Qew/meq];

  A = [A11 A12; A21 A22];
  B = [B1; B2];
  C=A(2,:);
  C(3)=0;
  C(4)=0;
  D=B(2);

  M = ss(A, B, C, D);
  G = tf(M);
  mag=bode(G,w);
  mag=squeeze(mag);
end