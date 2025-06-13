% This code aims to duplicate Figs. 2 and 3 in
%
%    Vanvincq and Bouwmans, "Numerical Simulation of the Nonlinear
%    Propagation of Ultrashort Pulses in Fiber Amplifiers Using
%    Time-Frequency Representations," J. Light. Technol. 41, 314-319 (2023)

close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1030e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
%sim.gpu_yes = false;
sim.save_period = 0.1;

% -------------------------------------------------------------------------

% Gain fiber
sim.progress_bar_name = 'Gain';
fiber.L0 = 0.8; % m; fiber length
fiber.MFD = 20; % um; mode-field diameter
fiber.n2 = 0;

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider the Raman
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE_propagate(fiber,sim); % load default parameters for "fiber" and "sim"

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 20; % um
gain_rate_eqn.cladding_diameter = 400; % um
gain_rate_eqn.core_NA = 0.065;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 976; % nm
gain_rate_eqn.absorption_to_get_N_total = 1.45; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 4.9; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^18; % the number of time points
time_window = 15e3; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( sim,gain_rate_eqn,{ifftshift(lambda,1),ifftshift(lambda,1)} );

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
rep_rate = 100; % Hz

tfwhm = 0.2; % ps
total_energy = 10e3; % nJ

pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

% Stretched method 1
stretched_duration = 500; % ps
%[~,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',stretched_duration,30*pi/180,sim.lambda0*1e9,t,prop_output.fields,1e-3/1000,100 );
% Stretched method 2
chirp_sign = 1; % for positive chirp (normal dispersive delay)
omega = ifftshift(2*pi*f,1);
spectrum_amplitude = ifft(prop_output.fields);
func_chirp = calc_chirp;
[chirp,stretched_field] = func_chirp.General( stretched_duration,omega,spectrum_amplitude,chirp_sign );

prop_output.fields = struct('forward', stretched_field,...
                            'backward',zeros(Nt,1));

prop_output_woPulse = struct('dt',(1/rep_rate*1e12-time_window)/Nt,...
                             'fields',struct('forward', zeros(Nt,1),...
                                             'backward',zeros(Nt,1)));

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output(2).fields.forward).^2,1),2)*prop_output(2).dt/10^3); % energy in nJ

t_woPulse = (-Nt/2:Nt/2-1)'*prop_output_woPulse.dt; % ps

save('one_nsPulse.mat');