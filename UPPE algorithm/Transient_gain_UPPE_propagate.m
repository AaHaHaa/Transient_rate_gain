function foutput = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn)
%TRANSIENT_GAIN_UPPE_PROPAGATE Propagate an initial pulse through an 
%arbitrary distance of an optical fiber with UPPE
%   This is a caller function, calling
%   UPPE_propagate_with_adaptive() or
%   UPPE_propagate_without_adaptive()
%   based on whether to use adaptive step-size method or not.
% -------------------------------------------------------------------------

% Curernt code supports only scalar and linearly-polarized simulations
sim.scalar = true; % non-polarized
sim.ellipticity = 0; % linearly polarized
sim.midx = 1; % single spatial mode

%%
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));

sim.cuda_dir_path = [upper_folder 'cuda'];

%% Determine whether to use adaptive-step-size method
% Don't use adaptive-step-size method if ASE is included or under counterpumping.
% Since back-and-forth iterations are needed for counterpumping or bidirectional pumping anyway, it'll just run with ASE.
% This line should be set in gain_info() already. Set it here again to be careful.
if gain_rate_eqn.counterpump_power ~= 0
    gain_rate_eqn.ignore_ASE = false;
    gain_rate_eqn.include_ASE = true;
end

% Single-pass adaptive-step simulations are fast.
if gain_rate_eqn.include_ASE
    adaptive_dz_str = 'without';
else
    adaptive_dz_str = 'with';
end

UPPE_propgation_func = str2func(['UPPE_propagate_', adaptive_dz_str, '_adaptive']);

%% Run the pulse propagation
foutput = UPPE_propgation_func(fiber, initial_condition, sim, gain_rate_eqn);

end