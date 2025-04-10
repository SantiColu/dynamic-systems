% condición de vuelo sin deslizamiento ni velocidad de cabeceo
%
% W  : peso
% y  : gamma [⁰]
% QS : presión dinámica x superficie de referencia [Pa]
% FLP: deflexión de flaps [⁰]
%
% trim
%     .cl
%     .cd
%     .alfa [⁰]
%     .alfa_idx 
%     .T  empuje [N]
%
function trim = level_flight(W, y, QS, data, FLP, flp_data)

    if nargin > 4 && FLP > 0
        flap_idx = interp1(flp_data.delta, 1:length(flp_data.dcl), FLP, 'nearest');
        cl_flap  = flp_data.dcl(flap_idx);
    else
        cl_flap  = 0;
    end

    % descartamos la zona de pérdida
    [~, k] = max(data.cl);
    data.cl    = data.cl(1:k);
    data.cd    = data.cd(1:k);
    data.alpha = data.alpha(1:k);

    % sustentación para equilibrio con la senda de planeo seleccionada
    Wy  = W * cosd(y); 
    % valores iniciales
    T   = Wy/20 + W * sind(y);  % empuje
    a   = 0; % ángulo de ataque [⁰]

    for k=1:100
        a_ = a;
        T_ = T;
        % coeficiente de sutentación requerido para esa condición de vuelo
        % suma de fuerzas en z wind = 0
        cl  = (Wy - T*sind(a)) / QS - cl_flap;
        % 'angulo de ataque para ese cl
        a   = (interp1(data.cl, data.alpha, cl) + a_)/2;
        if isnan(a)
            error 'vuelo nivelado no posible con esta combinación de carga y velocidad' 
        end    
        cd  = interp1(data.alpha, data.cd , a);    
        if cl_flap > 0
            alfa_idx = floor(interp1(data.alpha, 1:length(data.alpha), a));
            cd_flap  = flp_data.dcdmin(flap_idx) + flp_data.dcdi(alfa_idx, flap_idx);
            cl = cl + cl_flap;
            cd = cd + cd_flap;
        end
        ca  = cd * cosd(a) - cl * sind(a);
        D   = QS * ca;
        p   = (y + a) * pi / 180;
        T   = (W*sin(p) + D + T_)/2;
        if abs(a - a_) < 1e-3 && abs(T - T_) < 1e-1 
            break;
        end
    end

    trim.cl       = cl;
    trim.cd       = cd;
    trim.alfa     = a;
    trim.T        = T;
end
