%% Argumentos
%
% coef: derivativas de la aeronave
% elev: derivativas del elevador
% prop
%    .Tmax     : empuje máximo (si es reactor)
%    .P_rd     : propeller radius (si usa hélice)
%    .d_cg     : distancia del eje de empuje al CG  
% ref: condición de vuelo [unidades mks]
%    .speed    : TAS 
%    .q        : presión dinámica 
%    .g        : aceleración de la gravedad
%    .CL        
%    .CD     
%    .aoa     : ángulo de ataque [⁰]
%    .mass     
%    .inertia  : tesnor de inercia
%    .pitch    : ángulo de cabeceo [⁰]
%
% inp: 1: { elev thr }, 2: { uw ww qw }, otro: { elev thr uw ww qw },
% alfa_state: si es true segundo estado en el modelo se convierte a ángulo de ataque
%
%% Retorno
% mtx : matrices de parámetros
%       .M    : matriz de inercia
%       .Aa   : coeficientes aerodinámicos
%       .Ai   : acoplamientos de inercia
%       .Be   : sensibilidad al elevador (comando normalizado [-1 +1])
%       .Bt   : sensibilidad al empuje (comando normalizado [0 1])
%
% MS : modelo de estados
%
function [MS, mtx, aprox] = modelo_longitudinal(coef, elev, prop, ref, inp, alfa_state)

    if nargin < 5 || isempty(inp)
        inp = 0;
    end
    if nargin < 6 || isempty(alfa_state)
        alfa_state = true;
    end

    c = coef.cbar;
    S = coef.sref;

    Uo     = ref.speed;
    Q      = ref.q;
    g      = ref.g;
    CLo    = ref.CL;
    CDo    = ref.CD;
    aoa    = ref.alfa;
    pitch  = ref.pitch;
    gamma  = pitch - aoa;
    m      = ref.mass;
    Iyy    = ref.inertia(2,2);
   
    if isfield(ref, 'xcg')
        x_ref = -ref.xcg; % posición del CG en coordenadas body
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


    %% COEFICIENTES AERODINAMICOS (wind axes)
    CDu = 0; 
    CLu = 0; 

    if ~isfield(coef, 'cda')
        coef.cda = derivada_parabolica(coef.alpha, coef.cd');% 1/rad
    end
    CDa    = interp1(coef.alpha, coef.cda , aoa, 'pchip')*fderiv; % 1/rad
    CLa    = interp1(coef.alpha, coef.cla , aoa, 'pchip')*fderiv; % 1/rad
    CLadot = interp1(coef.alpha, coef.clad, aoa, 'pchip')*fderiv; % 1/rad
    
    Cmu    = 0; 
    Cmadot = interp1(coef.alpha, coef.cmad, aoa, 'pchip')*fderiv; % 1/rad
 
    k = find(isnan(coef.cma),1);
    if ~isempty(k)
        coef.cma(k:end) = coef.cma(k-1);
    end
    Cma    = interp1(coef.alpha, coef.cma , aoa, 'pchip')*fderiv; % 1/rad
    Cmq    = interp1(coef.alpha, coef.cmq , aoa, 'pchip')*fderiv; % 1/rad

    % Elevador
    if ~isempty(elev)
        dmax  = max(abs(elev.delta)); 
        delta = elev.delta/dmax;
        dcdi  = interp2(delta, coef.alpha, elev.dcdi, delta, aoa);
        cdde  = derivada_parabolica(delta, dcdi); % 1/rad
        clde  = derivada_parabolica(delta, elev.dcl'); % 1/rad
        cmde  = derivada_parabolica(delta, elev.dcm'); % 1/rad

        CDde = interp1(delta, cdde, 0, 'pchip');
        CLde = interp1(delta, clde, 0, 'pchip');
        Cmde = interp1(delta, cmde, 0, 'pchip');
        
        %
        % idx  = floor(interp1(elev.delta, 1:length(elev.delta), 0));
        % dlt  = (elev.delta(idx+1)-elev.delta(idx-1)) / dmax; % normalizado;
        % CDde = (elev.dcdi(a_idx, idx+1) - elev.dcdi(a_idx, idx-1))/dlt
        % CLde = (elev.dcl(idx+1) - elev.dcl(idx-1))/dlt
        % Cmde = (elev.dcm(idx+1) - elev.dcm(idx-1))/dlt
    end

    %% Derivativas para X, Z (wind axes) y M 
    QS    = Q*S;    % factor dimensional para coeficientes de fuerza 
    QSc   = QS*c;   % factor dimensional para coeficientes de momento
    W     = m*g;    % peso 
    w_    = 2*Uo/c; % velocidad angular de referencia
    
    Xu    = -(CDu+2*CDo/Uo)*QS;    % dF/dU = 2.Uo.F
    Xw    =  (CLo-CDa)/Uo  *QS;    % dF/dW = dF/da x da/dW = dF/da x 1/Uo
    Xwdot =   0;

    Zu    = -(CLu+2*CLo/Uo)*QS;
    Zw    = -(CLa+CDo)  *QS/Uo;
    Zq    =   0;
    Zwdot = -(CLadot/Uo)*QS/w_;

    Mu    = Cmu       * QSc;
    Mw    =(Cma/Uo)   * QSc; 
    Mq    = Cmq       * QSc/w_; 
    Mwdot =(Cmadot/Uo)* QSc/w_; 

    %% Derivativas para los comandos
    % Elevador (wind axes)
    if ~isempty(elev)
        Xde =-CDde*QS ;  
        Zde =-CLde*QS ;  
        Mde = Cmde*QSc;  
    else
        Xde = 0;
        Zde = 0; 
        Mde = 1;  
    end
    % Empuje (wind axes)
    if ~isempty(prop)
        p_pitch = aoa; % respecto de la dirección de vuelo
        if isfield(prop, 'pitch')
            p_pitch = p_pitch + prop.pitch;
        end
        St = sind(p_pitch);
        Ct = cosd(p_pitch); 
        Tu = QS*CDo + W*sind(gamma);
       %Tw = Tu*tan(p_pitch); 
        if isfield(prop, 'Tmax')
            Xdt =  Ct        * prop.Tmax;
            Zdt =  St        * prop.Tmax;
            Mdt = -prop.d_cg * prop.Tmax;   
        else
            To = Tu/cosd(p_pitch);
            %Po = To*Uo;
            R  = prop.P_rd;
            ef = 2 / (1 + sqrt(To/(Q*pi*R^2)+1));
            ef = 0.8 * ef;
            
            Tmax = prop.Pmax/Uo * ef;
            Xdt =  Ct        * Tmax;
            Zdt =  St        * Tmax;
            Mdt = -prop.d_cg * Tmax;   
        end
    else
        Xdt = 0;
        Zdt = 0;
        Mdt = 0;
        Tu  = 0;
       %Tw  = 0;
    end

    % Ajustes por corrimiento de del CG (a revisar)   
    Mu = Mu - x_ref * Zu; 
    Mw = Mw - x_ref * Zw;   
    % La velocidad de cabeceo produce un cambio en w = x_ref.q en el punto 
    % de referencia, y con ello la fuerza y momento asociado 
    Mq = Mq - x_ref * Mw; 
    Zq = Zq - x_ref * Zw;

    %% MODELO EN ESPACIO DE ESTADOS
    %  Como tenemos derivativas en función del ángulo de ataque, inicialmente
    %  planteamos todo en wind axes
    
    So = sind(gamma);
    Co = cosd(gamma);
    Aa = [Xu Xw 0  0
          Zu Zw Zq 0
          Mu Mw Mq 0 
          0  0  0  0];   
    Ai = [0   0   0    -W*Co 
          0   0   m*Uo -W*So %-Tu
          0   0   0     0 
          0   0   1     0];
    Be = [Xde ; Zde ; Mde ; 0 ];
    Bt = [Xdt ; Zdt ; Mdt ; 0 ];

    % INERCIA GENERALIZADA
    M  = [m  -Xwdot  0    0 
          0 m-Zwdot  0    0 
          0  -Mwdot  Iyy  0
          0     0    0    1];

    % verificaciones de período largo
    %
    aprox.lp.w =  sqrt(-Zu*g/m/Uo);
    aprox.lp.d =  CDo/(sqrt(2)*CLo); % = -Xu/(m*2*wlp)  
    aprox.sp.w =  sqrt(Zw*Mq/(m*Iyy) - Mw*Uo/Iyy);
    aprox.sp.d = -((Mq+Mwdot*Uo)/Iyy+(Zw)/m)/(2*aprox.sp.w);
    
    A_ = M\(Aa+Ai);
      
    %% De wind axes a body axes
    %
    % |xw| = | Ca  Sa| |xb|
    % |zw|   |-Sa  Ca| |zb|
    %
    % cd =  ca.Ca + cn.Sa
    % cl = -ca.Sa + cn.Ca
    %
    % Verificación:
    %
    % clf; a = coef.alpha'; sn = sind(a); cs = cosd(a); 
    % subplot(211); plot(a, coef.cd); hold on; grid on; plot(a, coef.ca.*cs+coef.cn.*sn); 
    % subplot(212); plot(a, coef.cl); hold on; grid on; plot(a,-coef.ca.*sn+coef.cn.*cs);
    %   
    Sa = sind(aoa);
    Ca = cosd(aoa);               
    Twb = [Ca  Sa  0   0  
          -Sa  Ca  0   0     
           0   0   1   0 
           0   0   0   1 ]; 
   
    mtx.M   = M *Twb;
    mtx.Aa  = Aa*Twb;
    mtx.Ai  = Ai*Twb;
    mtx.Be  = Be;
    mtx.Bt  = Bt;
    mtx.Twb = Twb;
    
    A = mtx.M\(mtx.Aa+mtx.Ai);
    
    % verificación del cambio de coordenadas
    assert(max(abs(sort(eig(A_))-sort(eig(A)))) < 1e-6);
    
    %%
    switch inp
        case 1
            if isempty(prop)
                entradas = {'elv'};
                B = mtx.Be; 
            else
                entradas = {'elv' 'thrust'};
                B = [mtx.Be mtx.Bt]; 
            end
        case 2
            entradas = {'uw' 'ww' 'qw'};
            B = mtx.Aa(:,1:3);
        otherwise
            if  isempty(prop)
                entradas = {'elv' 'uw' 'ww' 'qw'};
                B = [mtx.Be mtx.Aa(:,1:3)];
            else
                entradas = {'elv' 'thrust' 'uw' 'ww' 'qw'};
                B = [mtx.Be mtx.Bt mtx.Aa(:,1:3)];
            end
    end    
    B = mtx.M\B;	 

    %% cambio de coordenas aoa en lugar de w
    if alfa_state
        T = diag([1 Uo 1 1]);
        A = T\A*T;
        B = T\B;
        estados  = {'u' 'alpha' 'q' 'pitch'}; 
    else
        estados  = {'u' 'w' 'q' 'pitch'}; 
    end
    salidas = estados;

%     C = eye(length(A));
%     D = zeros(length(A),size(B,2));
% %    if exist("OCTAVE_VERSION", "builtin") 
% %      MS = ss(A, B); %, 'STNAME' , estados, 'INNAME' , entradas, 'OUTNAME', salidas);
% %    else
%       MS = ss( A, B, C, D ...
%              , 'statename' ,estados ...
%              , 'inputname' ,entradas ...
%              , 'outputname',salidas);
% %    end      
    MS.A = A;
    MS.B = B;
    MS.StateName = estados;
    MS.InputName = entradas;

end