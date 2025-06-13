% This code demonstrates the nonlinear CPA.
%
% For details, please refer to
% Kuznetsova et al., "Chirped-pulse amplification near the gain-narrowing
% limit of Yb-doped fiber using a reflection grism compressor," Appl. Phys.
% B 88, 515-518 (2007).
%
% Other recommended papers of nonlinear CPA, please read
% 1. Zhou et al., "Compensation of nonlinear phase shifts with third-order
%    dispersion in short-pulse fiber amplifiers," Opt. Express 13, 4869-
%    4877 (2005).
% 2. Kuznetsova and Wise, "Scaling of femtosecond Yb-doped fiber amplifiers
%    to tens of microjoule pulse energy via nonlinear chirped pulse 
%    amplification," Opt. Lett. 32, 2671-2673 (2007).

close all; clearvars;

addpath('../../../../../../GMMNLSE/GMMNLSE algorithm/','../../../../../../GMMNLSE/user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_GMMNLSE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1030e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain';
fiber.L0 = 0.8; % m; fiber length
fiber.MFD = 10; % um; mode-field diameter
fiber.n2 = 0;
sim.gpuDevice.Index = 2;

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider the Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode'); % load default parameters for "fiber" and "sim"

sim.save_period = fiber.L0/10;

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% Activating "reuse_data" or "linear_oscillator_model" requires setting other parameters.
% Check the example or "gain_info.m".
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.08;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 1.7; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0.345; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/100; % Assume 100 Hz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^19; % the number of time points
time_window = 15e3; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% calculate fiber betas from silica refractive index
% This is important to correctly simulate the broadband situations.
% Taylor-series coefficients is only good in narrowband situations.

% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Setup initial conditions
pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
tfwhm = 0.08; % ps
total_energy = 15e3; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

stretched_duration = 1500; % ps
[~,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',stretched_duration,30*pi/180,sim.lambda0*1e9,t,initial_pulse.fields,1e-3/1000,100 );
initial_pulse.fields = stretched_field;

%% Run the simulations
prop_output = GMMNLSE_propagate(fiber, initial_pulse, sim, gain_rate_eqn); 

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output.fields).^2,1),2)*prop_output.dt/10^3); % energy in nJ

figure;
plot(prop_output.z,energy,'linewidth',2);
xlabel('Propagation distance (m)');
ylabel('Pulse energy (nJ)');
set(gca,'fontsize',20);

figure;
plot(t,abs(prop_output.fields(:,:,end)).^2/1e3,'linewidth',2,'Color','b');
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);

save('Steady-state.mat');