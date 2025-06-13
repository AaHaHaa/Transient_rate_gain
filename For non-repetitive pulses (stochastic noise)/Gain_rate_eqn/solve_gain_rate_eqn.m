function N = solve_gain_rate_eqn(sim,gain_rate_eqn,...
                                 N,...
                                 R,...
                                 Nt,dt)
%SOLVE_GAIN_RATE_EQN Solves the populations
%
% computational dimension: (Nt,num_cross_sections/num_levels)

% Load required helper functions
func = solve_gain_rate_eqn_helpers();

% -------------------------------------------------------------------------
% --------------------- Rate equation to get N ----------------------------
% -------------------------------------------------------------------------

% Ion density among various levels
N = func.solve_N(sim,gain_rate_eqn,...
                 gain_rate_eqn.acyclic_conv_stretch(Nt),dt,...
                 N,gain_rate_eqn.N_total,...
                 R); % unit: 1/um^3

% Downsample it back to the original time window
N = N(1:Nt,:);

end