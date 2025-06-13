function foutput = Periodic_transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn)
%PERIODIC_TRANSIENT_GAIN_UPPE_PROPAGATE Propagate an initial pulse through an 
%arbitrary distance of an optical fiber with UPPE incorporating transient
%gain
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

addpath([upper_folder,'Gain_rate_eqn']);

sim.cuda_dir_path = [upper_folder,'cuda'];

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

%% Apply the narrowband transformation (due to the scaled Fourier transform)
scaledFT_func = narrowband_scaledFT();

if sim.cs > 1
    num_windows = length(initial_condition);
    for window_i = 1:num_windows
        initial_condition(window_i).dt = initial_condition(window_i).dt*sim.cs;
        initial_condition(window_i).fields.forward = scaledFT_func.convert(initial_condition(window_i).fields.forward(:,:,end),sim.cs);
        if isfield(initial_condition(window_i).fields,'backward')
            initial_condition(window_i).fields.backward = scaledFT_func.convert(initial_condition(window_i).fields.backward(:,:,end),sim.cs);
        end
    end
end

%% Run the pulse propagation
foutput = UPPE_propgation_func(fiber, initial_condition, sim, gain_rate_eqn);

%% Recover from the narrowband transformation
if sim.cs > 1
    for window_i = 1:num_windows
        foutput(window_i).dt = foutput(window_i).dt/sim.cs;
        if ~isempty(foutput(window_i).fields.forward) % In simulations without ASE, its window contains no field, which we skip here.
            foutput(window_i).fields.forward = scaledFT_recover_field(scaledFT_func,foutput(window_i).fields.forward,sim.cs);
        end
        if isfield(foutput(window_i).fields,'backward') % Only simulations with ASE has backward fields.
            if ~isempty(foutput(window_i).fields.backward) % Just a careful check, but this line should be removable.
                foutput(window_i).fields.backward = scaledFT_recover_field(scaledFT_func,foutput(window_i).fields.backward,sim.cs);
            end
        end

        % Recover populations
        foutput(window_i).population = scaledFT_recover_population(foutput(window_i).population,sim.cs);
    end
end

end

%% Helper functions the narrowband transformation due to the scaled Fourier transform
function B = scaledFT_recover_field(scaledFT_func,A,cs)

Nt = size(A,1);
Nz = size(A,3);

B = zeros(Nt*cs,1,Nz);
for zi = 1:Nz
    B(:,:,zi) = scaledFT_func.recover(A(:,:,zi),cs);
end

end

function rN = scaledFT_recover_population(N,cs)

Nt = size(N,1);
num_levels = size(N,2);
Nz = size(N,3);

recover_idx = linspace(1,Nt+1,Nt*cs+1)';
recover_idx = recover_idx(1:end-1); % remove the last "Nt+1" point

rN = zeros(Nt*cs,num_levels,Nz);
for zi = 1:Nz
    for lidx = 1:num_levels
        rN(:,lidx,zi) = interp1((1:Nt)',N(:,lidx,zi),recover_idx,'linear','extrap');
    end
end

end