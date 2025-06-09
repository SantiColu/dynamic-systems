function [M,K] = femGeneral(pl)
  if nargin < 1
    pl = false; % Por defecto no plotear
  end

  [Mt, Kt, xt] = femTorre();    % Matrices torre
  [Mp, Kp, rp] = femPala();     % Matrices pala
  masaNacela = 24000; % Masa de la nacela

  Npalas = 3;
  Ndof_nodo = 2;                % 2 DOFs por nodo (viga: desplazamiento y rotación)

  Lt = size(Mt, 1);             % DOFs torre
  Lp = size(Mp, 1);             % DOFs una pala

  % Cada pala comparte su primer nodo (2 DOFs) con el último nodo de la torre
  % Entonces cada pala contribuye con Lp - 2 DOFs únicos
  DOFs_extra = Npalas * (Lp - Ndof_nodo);
  TotalSize = Lt + DOFs_extra;

  % Inicializamos matrices globales
  M = zeros(TotalSize);
  K = zeros(TotalSize);

  % Agrego la masa de la nacela
  % TODO: REVISAR POR QUE ME BAJA TANTO EL PRIMER MODO (ESTA OK?)
  % Mt(end-1:end-1) = Mt(end-1:end-1) + masaNacela;

  % Ensamblar torre
  M(1:Lt, 1:Lt) = Mt;
  K(1:Lt, 1:Lt) = Kt;

  if pl
    plotMatrix(1, K);
  end

  % Ensamblar palas
  for i = 0:Npalas-1
    % DOFs locales de la pala
    dofs_cabeza = 1:Ndof_nodo;                        % Primer nodo (conectado a torre)
    dofs_cuerpo = Ndof_nodo+1:Lp;                     % Resto de la pala

    % DOFs globales
    dofs_torre = Lt - Ndof_nodo + 1 : Lt;             % Último nodo de la torre
    st = Lt + i * (Lp - Ndof_nodo) + 1;               % Inicio bloque de pala en global
    dofs_pala_global = st : st + length(dofs_cuerpo) - 1;

    % Parte propia de la pala
    M(dofs_pala_global, dofs_pala_global) = Mp(dofs_cuerpo, dofs_cuerpo);
    K(dofs_pala_global, dofs_pala_global) = Kp(dofs_cuerpo, dofs_cuerpo);

    % Nodo compartido (acople torre ↔ pala)
    M(dofs_torre, dofs_torre) += Mp(dofs_cabeza, dofs_cabeza);
    K(dofs_torre, dofs_torre) += Kp(dofs_cabeza, dofs_cabeza);

    M(dofs_torre, dofs_pala_global) += Mp(dofs_cabeza, dofs_cuerpo);
    M(dofs_pala_global, dofs_torre) += Mp(dofs_cuerpo, dofs_cabeza);

    K(dofs_torre, dofs_pala_global) += Kp(dofs_cabeza, dofs_cuerpo);
    K(dofs_pala_global, dofs_torre) += Kp(dofs_cuerpo, dofs_cabeza);

    if pl 
      plotMatrix(i+2, K);
    end
  end

  % Empotrar base torre (primer nodo)
  M = M(Ndof_nodo+1:end, Ndof_nodo+1:end);
  K = K(Ndof_nodo+1:end, Ndof_nodo+1:end);

  % removemos los giros en el nodo de union
  % TODO: REVISAR POR QUE ME SACA EL PRIMER MODO (ESTA OK?)
  % K(dofs_torre,:) = [];
  % K(:, dofs_torre) = [];
  % M(dofs_torre,:) = [];
  % M(:, dofs_torre) = [];
end 
