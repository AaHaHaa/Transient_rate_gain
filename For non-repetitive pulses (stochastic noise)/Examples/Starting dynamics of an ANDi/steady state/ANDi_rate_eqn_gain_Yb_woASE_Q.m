% This code finds the mode-locked state of an all-normal dispersion (ANDi) 
% oscillator.
% An ANDi laser typically contains a sequence of a long passive fiber, a 
% short gain fiber, and a short passive fiber, following by an output
% coupler, a saturable absorber, and a spectral filter.
%
% This file uses the rate-equation gain model for the gain fiber.

close all; clearvars;

addpath('../../../GMMNLSE/GMMNLSE algorithm/','../../../GMMNLSE/user_helpers/');
    
%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore , the backward-propagating pulses need to be taken into account.
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 6; % um
gain_rate_eqn.core_NA = 0.2;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 280; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0.3; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % For counterpumping or considering ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% General parameters
sim.lambda0 = 1030e-9; % m; center wavelength
sim.f0 = 2.99792458e-4/sim.lambda0; % THz; center frequency
%sim.progress_bar = false;
sim.gpu_yes = false; % use GPU or not

sim6 = sim; % for the passive fibers
sim6.progress_bar_name = 'SMF (6um)';
sim6.save_period = 1e-3; % m

[fiber6,sim6] = load_default_GMMNLSE_propagate([],sim6,'single_mode'); % for the passive fibers

% -------------------------------------------------------------------------
% Long SMF
fiber_6a = fiber6;
fiber_6a.L0 = 10;

% Gain fiber
sim_Gain_6 = sim;
sim_Gain_6.gain_model = 2; % use the rate-eqn-gain model
sim_Gain_6.progress_bar_name = 'Gain (6um)';
sim_Gain_6.save_period = 1e-3; % m
fiber_Gain_6.L0 = 0.2;
[fiber_Gain_6,sim_Gain_6] = load_default_GMMNLSE_propagate(fiber_Gain_6,sim_Gain_6,'single_mode'); % for the gain fibers
fiber_Gain_6.betas = fiber_Gain_6.betas*0;
fiber_Gain_6.n2 = 0;

% Short SMF
fiber_6b = fiber6;
fiber_6b.L0 = 1;

% -------------------------------------------------------------------------
% ----------------------------------- All ---------------------------------
% -------------------------------------------------------------------------
fiber_cavity = [fiber_Gain_6];
sim_cavity = [sim_Gain_6];

%% Setup general cavity parameters
max_rt = 500; % maximum roundtrips (in case it doesn't converge)
Nt = 2^11; % the number of points
time_window = 1e5; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
OC = 0.8; % output coupling
loss = 0.5; % the total loss of the cavity
saturation_power = 250; % the saturation power of the saturable absorber; W
moddepth = 0.9; % the modulation depth of the saturable absorber
tol_convergence = 1e-3; % the tolerance of the convergence of the output ANDi pulse energy

%% Spectral filter parameters
gaussexpo = 1;
plot_filter_yes = true;
spectral_filter = struct('bw',6, ... % bandwidth (nm)
                         'cw',1028); % center wavelength (nm)

%% Setup initial conditions
tfwhm = 1e4; % ps
total_energy = 1; % nJ

prop_output = build_MMgaussian(tfwhm, time_window, total_energy,1,Nt);

%% Saved field information
L0 = sum([fiber_cavity.L0]); % total fiber length
save_num = sum(int64([fiber_cavity.L0]./[sim_cavity.save_period])) + 1; % the total number of fields to save
save_num = double(save_num);
saved_z = zeros(1,save_num); % the propagation distance to save
field = cell(1,max_rt);
pump = cell(1,max_rt); % the pump power to save
splice_z = cumsum([fiber_cavity.L0]); % where the splice points are
output_field = zeros(Nt,1,max_rt); % the output pulse
max_save_per_fiber = 20;

%% Load gain parameters
L_air = 29.69; % 1 is the free-space length
c = 299792458; % m/s
v = 1/fiber_cavity(1).betas(2)*1e12; % velocity in the fiber

t_rep = L0/v + L_air/c; % s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                        % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.t_rep = t_rep;

gain_rate_eqn = gain_info( fiber_Gain_6,sim_Gain_6,gain_rate_eqn,ifftshift(lambda,1) );

%% Run the cavity simulation

func = analyze_sim; % the container of several analyzing functions

% Initialize some parameters
output_energy = zeros(max_rt,1);
rt_num = 0;
pulse_survives = true;
while rt_num < max_rt
    time_delay = 0;
    current_z = 0;
    zn = 1;
    rt_num = rt_num + 1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);
    % -----------------------------------------------------------------
    for j = 1
        prop_output = GMMNLSE_propagate(fiber_cavity(j), prop_output, sim_cavity(j),gain_rate_eqn);
        prop_output_backup = prop_output;
            
        time_delay = time_delay + prop_output.t_delay(end);
            
        % Save the information
        [saved_field,saved_z_this_fiber] = func.extract_saved_field(prop_output.fields,max_save_per_fiber,current_z,prop_output.z);
        field{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_field;
        saved_z(zn:zn+size(saved_field,3)-1) = saved_z_this_fiber;
        
        current_z = saved_z_this_fiber(end);
        zn = zn + size(saved_field,3)-1;
    end

    saved_z = saved_z(1:zn);
    
    % -----------------------------------------------------------------
    % Output coupler
    output_field(:,:,rt_num) = sqrt(OC)*prop_output.fields(:,:,end);
    prop_output.fields = sqrt(1-OC)*prop_output.fields(:,:,end);
    
    % -----------------------------------------------------------------
    % Saturable absorber
    %a=prop_output;
    %prop_output = saturable_absorber_action_simple(prop_output, saturation_power, moddepth);
    
    % -----------------------------------------------------------------
    % Spectral filter
    if rt_num ~= 1
        %close(fig_filter); % close the previous figure and plot a new one
    end
    %[prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter
    
    prop_output.fields = sqrt(1-loss)*prop_output.fields(:,:,end);
    
    % -----------------------------------------------------------------
    % Energy of the output field
    output_energy(rt_num) = trapz(abs(output_field(:,:,rt_num)).^2)*prop_output.dt/1e3; % energy in nJ

    % If the energies stop changing, then we're done!
    if rt_num ~= 1
        close(fig); % close the previous figure and plot a new one
    end
    warning('off')
    output_field(:,:,rt_num) = pulse_tracker(output_field(:,:,rt_num));
    warning('on');
    [converged_yes,fig] = check_convergence( output_energy,output_field(:,:,rt_num),f,t,tol_convergence,true );
    
    % ---------------------------------------------------------------------
    % Display running time of each roundtrip
    t_iteration_end = toc(t_iteration_start);
    t_iteration_spent = datevec(t_iteration_end/3600/24);
    fprintf(': %1u:%2u:%3.1f\n',t_iteration_spent(4),t_iteration_spent(5),t_iteration_spent(6));
    
    % ---------------------------------------------------------------------
    % Update the repetition rate based on "time_delay"
    % The pulse will shift in time with respect to the moving frame.
    % Since I implement the pulse centering function in the code,
    % I can use the information of "time_delay" in each roundtrip to calibrate to get the actual repetition rate.
    gain_rate_eqn.t_rep = t_rep + time_delay*1e-12;
    
    % ---------------------------------------------------------------------
    % Characteristic lengths
    [ dispersion_length,nonlinear_length ] = characteristic_lengths( abs(output_field(:,1,rt_num)).^2,t,sim.f0,fiber_cavity(end).betas(3,:),1/fiber_cavity(end).SR );
    fprintf('  L_DL=%4f(m)\n',dispersion_length);
    fprintf('  L_NL=%4f(m)\n',nonlinear_length);
    
    % ---------------------------------------------------------------------
    % Plot
    if rt_num ~= 1
        close(fig_evolution); % close the previous figure and plot a new one
    end
    fig_evolution = func.analyze_fields(t,f,field{rt_num},saved_z,splice_z);
    
    % ---------------------------------------------------------------------
    % Break if converged
    if converged_yes
        cprintf('blue','The field has converged!\n');
        break;
    end
    % Break if pulse dies
    if output_energy(rt_num) < 0.01 % smaller than 0.01 nJ
        disp('The pulse dies.');
        pulse_survives = false;
        break;
    end
end

%% Finish the simulation and save the data
% Clear reducdant parts of the data
field = field(1:rt_num);
output_field = output_field(:,:,1:rt_num);
energy = output_energy(arrayfun(@any,output_energy)); % clear zero

%close(fig,fig_filter,fig_evolution);
close all;

%% Compress the pulse
[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,output_field(:,:,end),'Treacy-t',pi/6,1e-3/1000,true );

save('ss_for_initial_N_Q.mat');