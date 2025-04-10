function mdl = A320_build_model(mass_frac, height, knt, flap, save_result, alfa_state, gam, slp)
format long
    if nargin < 1 || isempty(mass_frac)
        mass_frac = 1;
    end
    if nargin < 2 || isempty(height)
        height = 1000;
    end
    if nargin < 3 || isempty(knt)
        speed = 500/3.6;
    else
        speed = knt*1.852/3.6;  % knt -> m/s 
    end
    if nargin < 4 || isempty(flap)
        flap = 0;
    end
    if nargin < 5 || isempty(save_result)
        save_result = false;
    end
    if nargin < 5
        alfa_state = false;
    end

    modelo = 'A320';
    load(strcat('./', modelo, '_data'));
    old_path = addpath('./core');
    
  
        mdl.gamma   = gam; % [⁰]
        % mdl.mass    = data.mass.zfw + (data.mass.max_takeoff - data.mass.zfw) * mass_frac;
        mdl.mass    = data.mass.max_landing;
        mdl.inertia = data.mass.inertia * mdl.mass / data.mass.max_takeoff;
        mdl.height  = height;
        mdl.speed   = speed; 
        mdl.flaps   = flap;
        mdl.xcg     = 0;
        mdl.rho     = densidad(mdl.height);     % kg/m3 
        mdl.g       = 9.81;                  % m/s²
        mdl.q       = 0.5*mdl.rho*mdl.speed^2;  % Presion dinamica [Pa]
    
        % mdt(1) = mdl.mass*mdl.g;
        % mdt(2) = mdl.gamma;
        % mdt(3) = mdl.q*data.sref;
        % % mdt(4) = data;
        % mdt(5) = mdl.flaps;
        % % mdt(6) = data.flap;
        % % data.flap

        % mdt
        
        trim = level_flight(mdl.mass*mdl.g, mdl.gamma, mdl.q*data.sref, data, mdl.flaps, data.flap);
    
        mdl.alfa     = trim.alfa;
        mdl.CL       = trim.cl;
        mdl.CD       = trim.cd; 
        mdl.pitch    = mdl.alfa + mdl.gamma;
    
        [mdl.lng, mdl.mtx.lng] = modelo_longitudinal(data, data.elev, data.thrust, mdl, [], alfa_state);
        [mdl.lat, mdl.mtx.lat] = modelo_lateral     (data, data.aler, data.rudr, mdl);
    
        path(old_path);
        if save_result
            model_file_name = strcat( modelo ...
                                    , '_g=', num2str(mdl.gamma) ...
                                    , '_h=', num2str(mdl.height) ...
                                    , '_v=', num2str(mdl.speed) ...
                                    , '_f=', num2str(mdl.flaps) ...
                                    , '.mat')
            save(model_file_name, 'mdl')
        end

    path(old_path);
end





