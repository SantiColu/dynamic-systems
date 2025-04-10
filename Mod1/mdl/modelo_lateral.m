% m: masa
% I: tensor de inercia (mks)
% H: altitud [m]
% Vo: TAS [m/s]
% FLP: posición del flap [°]
% coef: datos del datcom 
% flp_data: datos del datcom con flaps
% b : envergadura [m]
% c : cuerda de referencia [m]
% S : superficie de referencia [m2]
function [MS, mtx] = modelo_lateral(coef, aler, rudr, ref, beta_state)

    if nargin < 5 || isempty(beta_state)
        beta_state = true;
    end

    b  = coef.blref;
    S  = coef.sref;
    
    Vo    = ref.speed;
    Q     = ref.q;
    g     = ref.g;
    alfa  = ref.alfa;

    We = Vo*sind(alfa);
    Ue = Vo*cosd(alfa);
    Co = cosd(ref.pitch);
    
    m   = ref.mass;
    Ixx = ref.inertia(1,1);
    Ixz = ref.inertia(1,3);
    Izz = ref.inertia(3,3);
    
    if isfield(ref, 'xcg') 
        x_ref = -ref.xcg;
    else
        x_ref = 0;
    end
    switch coef.deriv
        case 'rad'
            fderiv = 1;
        case 'deg'
            fderiv = 180/pi;
        otherwise
            fderiv = 1;
    end    

    %% COEFICIENTES AERODINAMICOS
    Cyv  = interp1(coef.alpha, coef.cyb, alfa)*fderiv/Vo;
    Cyp  = interp1(coef.alpha, coef.cyp, alfa)*fderiv;
    Cyr  = 0; % no hay
    Clv  = interp1(coef.alpha, coef.clb, alfa)*fderiv/Vo;
    Clp  = interp1(coef.alpha, coef.clp, alfa)*fderiv;
    Clr  = interp1(coef.alpha, coef.clr, alfa)*fderiv;
    Cnv  = interp1(coef.alpha, coef.cnb, alfa)*fderiv/Vo;
    Cnp  = interp1(coef.alpha, coef.cnp, alfa)*fderiv;
    Cnr  = interp1(coef.alpha, coef.cnr, alfa)*fderiv;
    % Alerón
    if ~isempty(aler) 

        dmax      = max(abs(aler.delta)); 
        delta     = aler.delta/dmax;
        cnda      = interp2(delta, coef.alpha, aler.cn, delta, alfa);

        cnda = derivada_parabolica(delta, cnda); 
        clda = derivada_parabolica(delta, aler.clroll'); 

        Cnda = interp1(delta, cnda, 0, 'pchip');
        Clda = interp1(delta, clda, 0, 'pchip');

        %     idx  = floor(interp1(aler.delta, 1:length(aler.delta), 0));
        %     dlt  = (aler.delta(idx+1)-aler.delta(idx-1)) / max(abs(aler.delta)); % normalizado
        %     Cnda = (aler.cn(a_idx, idx+1) - aler.cn(a_idx, idx-1))/dlt;
        %     Clda = (aler.clroll(idx+1) - aler.clroll(idx-1))/dlt;
    end
    % Rudder
    if ~isempty(rudr)    
        dmax      = max(abs(rudr.delta)); 
        delta     = rudr.delta/dmax;

        cndr = derivada_parabolica(delta, rudr.dcn); 
        cydr = derivada_parabolica(delta, rudr.dcy); 
        cldr = derivada_parabolica(delta, rudr.clroll); 
        
        Cydr = interp1(delta, cydr, 0, 'pchip');
        Cldr = interp1(delta, cldr, 0, 'pchip');
        Cndr = interp1(delta, cndr, 0, 'pchip');

        %     idx  = floor(interp1(rudr.delta, 1:length(rudr.delta), 0));
        %     dlt  = (rudr.delta(idx+1)-rudr.delta(idx-1)) / max(abs(rudr.delta)); % normalizado
        %     Cydr = (rudr.dcy(idx+1) - rudr.dcy(idx-1))/dlt;
        %     Cndr = (rudr.dcn(idx+1) - rudr.dcn(idx-1))/dlt;
    end

    %% DEFINICION DE DERIVATIVAS A UTILIZAR
    QS   = Q  *S;
    QSb  = QS *b;
    QSb2Vo = QSb*b/(2*Vo);
    % Se desprecia Yda
    Yv  = Cyv*QS;
    Lv  = Clv*QSb;
    Nv  = Cnv*QSb;
    Yp  = Cyp*QSb2Vo/b ;
    Lp  = Clp*QSb2Vo;
    Np  = Cnp*QSb2Vo;
    Yr  = Cyr*QSb2Vo;
    Lr  = Clr*QSb2Vo;
    Nr  = Cnr*QSb2Vo;
    % Alerón
    if ~isempty(aler)    
        Lda = Clda*QSb;
        Nda = Cnda*QSb;
    else    
        Lda = 1;
        Nda = 0;
    end
    if ~isempty(rudr)    
        % Rudder
        Ydr = Cydr*QS;
        Ldr = Cldr*QSb;
        Ndr = Cndr*QSb;
    else    
        Ydr = 0;
        Ldr = 0;
        Ndr = 1;
    end
    mtx.Ba = [0   ; Lda ; Nda ; 0];
    %% MODELO EN ESPACIO DE ESTADOS
    % LATERAL-DIRECCIONAL
    mtx.M = [m    0    0    0
             0    Ixx -Ixz  0 
             0   -Ixz  Izz  0 
             0    0    0    1];
    
    Nv = Nv + x_ref * Yv;  
    Np = Np + x_ref * Yp;   
    Nr = Nr + x_ref * Yr;   
     
    mtx.Aa= [ Yv   Yp   Yr   0
              Lv   Lp   Lr   0 
              Nv   Np   Nr   0
              0    0    0    0];
    mtx.Ai= [ 0    m*We -m*Ue  m*g*Co
              0    0     0     0 
              0    0     0     0
              0    1     0     0];
    mtx.Ba = [0   ; Lda ; Nda ; 0]; 
    mtx.Bw = mtx.Aa(:,1:3); % v, p, r 
    if isempty(rudr) 
        B = [mtx.Ba mtx.Bw];
        entradas = {'ail' 'vw'};
    else
        mtx.Br = [Ydr ; Ldr ; Ndr ; 0];
        B = [mtx.Ba mtx.Br mtx.Bw];
        entradas = {'ail' 'rdr' 'vw' 'pw' 'rw'};
    end
    A = mtx.M\(mtx.Aa+mtx.Ai);  
    B = mtx.M\B;

    %% cambio de coordenas beta en lugar de v
    if beta_state
        T = diag([Vo 1 1 1]);
        A = T\A*T;
        B = T\B;
        estados  = {'beta' 'p' 'r' 'roll'};
    else
        estados  = {'v' 'p' 'r' 'roll'};  
    end
    % C = eye(4);
    % D = zeros(4,size(B,2));
    % salidas = estados;
    % MS = ss( A, B, C, D ...
    %        , 'statename' , estados  ...
    %        , 'inputname' , entradas ...
    %        , 'outputname', salidas);
    MS.A = A;
    MS.B = B;
    MS.StateName = estados;
    MS.InputName = entradas;
end

