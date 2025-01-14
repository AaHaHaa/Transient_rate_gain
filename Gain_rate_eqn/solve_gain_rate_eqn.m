function N = solve_gain_rate_eqn(sim,gain_rate_eqn,...
                                 N,...
                                 Aw_forward,Aw_backward,...
                                 Power_pump_forward,Power_pump_backward,...
                                 cross_sections,E_photon,...
                                 Nt,dt)
%SOLVE_GAIN_RATE_EQN Solves the populations
%
% computational dimension: (Nt,num_cross_sections/num_levels)

% Load required helper functions
func = solve_gain_rate_eqn_helpers();

% Compute F^(-1)[A(omega) * sqrt(sigma/hbar/omega)]
At_forward = fft(Aw_forward.*sqrt(cross_sections./E_photon),[],1); % sqrt(m^2/s); size: (Nt,num_cross_sections)
% Compute F^(-1)[A(omega) * sqrt(sigma/hbar/omega)]
At_backward = fft(Aw_backward.*sqrt(cross_sections./E_photon),[],1); % sqrt(m^2/s); size: (Nt,num_cross_sections)

A2 = abs(At_forward).^2 + abs(At_backward).^2; % A2(t); the number of scattered photons

% -------------------------------------------------------------------------
% --------------------- Rate equation to get N ----------------------------
% -------------------------------------------------------------------------

% Ion density among various levels
N = func.solve_N(sim,gain_rate_eqn,...
                 gain_rate_eqn.acyclic_conv_stretch(Nt),dt,...
                 N,gain_rate_eqn.N_total,...
                 gain_rate_eqn.overlap_factor,gain_rate_eqn.cross_sections_pump,...
                 A2,...
                 Power_pump_forward,Power_pump_backward); % unit: 1/um^3

% Downsample it back to the original time window
N = N(1:Nt,:);

end