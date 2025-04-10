function rho = densidad(H)
    kapha = -0.000022558; %/m
    rho_sl=1.225; % a nivel del mar kg/m3
    rho = rho_sl*((1+kapha*H)^(4.2561)); % a la altura dada kg/m3
end