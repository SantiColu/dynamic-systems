function [t, u, f, S] = getTurbulence(V0)
  U0=V0;
  s=10;
  L=500;

  f=0:0.005:4;
  w=2*pi*f;

  S=((s^2*L)/(pi*U0)*((1+3*L^2*(w/U0).^2)./((1+L^2*(w/U0).^2).^2)))';

  k=s*sqrt(3*U0/(pi*L));
  z=U0/sqrt(3*L);
  p=U0/L;
  G=tf(k*[1 z],[1 2*p p^2]);

  dt=0.1;
  T=60;
  t=0:dt:T;

  rng(2);
  o=randn(size(t));
  u=lsim(G,o,t);
end