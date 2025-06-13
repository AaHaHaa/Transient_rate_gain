function N = solve_gain_rate_eqn_steadyState(gain_rate_eqn,gpu_yes,...
                                             At_forward,At_backward,...
                                             pump_power,...
                                             Nt,dt)
%SOLVE_GAIN_RATE_EQN_STEADYSTATE Solves the steady-state populations
%
% computational dimension: (Nt,num_cross_sections/num_levels,save_points,num_windows)

% Load required helper functions
func = solve_gain_rate_eqn_helpers();

if gpu_yes
    Aw_forward = ifft(cell2mat(cellfun(@gather,At_forward(:,:,:,2:2:end),'UniformOutput',false)),[],1);
else
    Aw_forward = ifft(cell2mat(At_forward(:,:,:,2:2:end)),[],1);
end
Aw_forward = Aw_forward*permute(sqrt((Nt*dt(2))/sum(Nt*dt)),[1,3,4,2]); % change it to the correct physical unit of W
                                                     % A2*time_window = pulse energy
                                                     % pulse energy/t_rep = W
                                                     % t_rep: the time for a pulse to finish one round trip
if ~isempty(At_backward)
    if gpu_yes
        Aw_backward = ifft(cell2mat(cellfun(@gather,At_backward(:,:,:,2:2:end),'UniformOutput',false)),[],1);
    else
        Aw_backward = ifft(cell2mat(At_backward(:,:,:,2:2:end)),[],1);
    end
    Aw_backward = Aw_backward*permute(sqrt((Nt*dt(2))/sum(Nt*dt)),[1,3,4,2]);
else
    Aw_backward = 0;
end

A2  = sum(abs(Aw_forward).^2 + abs(Aw_backward).^2,4);

% -------------------------------------------------------------------------
% --------------------- Rate equation to get N ----------------------------
% -------------------------------------------------------------------------

% Ion density in various levels
N = func.solve_N_SS(gain_rate_eqn,gpu_yes,...
                    gain_rate_eqn.N_total,...
                    gain_rate_eqn.E_photon{2},...
                    gain_rate_eqn.overlap_factor,gain_rate_eqn.cross_sections_pump,gain_rate_eqn.cross_sections{2},...
                    pump_power,...
                    A2); % unit: 1/um^3

end