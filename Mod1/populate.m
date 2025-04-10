
h_rg   = 5000:1000:5000;
v_rg = 130:10:140;
d_rg = -3:1:5;
slp_rg = -1:1:3;
gam_rg = -1:1:3;

length(h_rg)*length(v_rg)*length(d_rg)*length(slp_rg)*length(gam_rg)

for h = h_rg
  for V = v_rg
    for d = d_rg
      for slp = slp_rg
        for gam = gam_rg

          simName = sprintf("output/flight_sim_h%.0f_v%.1f_gam%.1f_slp%.1f_d%.1f.mat", h, V, gam, slp, d)
          if exist(simName, 'file')
              continue
          else
            [t, x] = run_flight_sim(h, V, gam, slp, d);
            save(simName, 'x', 't');
          end

        endfor
      endfor
    endfor
  endfor
endfor
