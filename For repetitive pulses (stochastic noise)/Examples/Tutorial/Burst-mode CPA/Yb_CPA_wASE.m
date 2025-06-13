% This code plots the schematic of transient rate-gain model.
%
% ASE is included:
%   (1) gain_rate_eqn.ignore_ASE = false
%   (2) sim.incoherent needs to be supplied for incoherent-field windows, 
%       which is then passed into gain_info to compute some parameters to save in gain_rate_eqn.

close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1030e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
%sim.gpu_yes = false;
sim.save_period = 0.01;
sim.dz = 1e-3; % m

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain';
fiber.L0 = 1; % m; fiber length
fiber.MFD = 6; % um; mode-field diameter

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
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.12;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.55; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 3; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 2; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^16; % the number of time points
pulse_time_window = 4e3; % ps
dt = pulse_time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

sim.incoherent.dt = dt/3;
incoherent_lambda0 = 1070e-9; % m
sim.incoherent.f0 = 2.99792458e-4/incoherent_lambda0;
incoherent_f = sim.incoherent.f0+(-Nt/2:Nt/2-1)'/(Nt*sim.incoherent.dt); % THz
incoherent_lambda = c./(incoherent_f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
sim.cs = 8; % apply narrowband transformation for the coherent field to speed up the simulation
gain_rate_eqn = gain_info( sim,gain_rate_eqn,{ifftshift(incoherent_lambda,1),ifftshift(lambda,1)} );

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
burst_rep_rate = 1e6; % Hz
intra_burst_rate = 100e6; % Hz; intra-burst rate

num_pulses = 10;

tfwhm = 0.4; % ps
total_energy = 0.1; % nJ

pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, pulse_time_window, total_energy, 1, Nt, {'ifft',freq_shift});

% Stretched method 1
stretched_duration = 500; % ps
[~,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',stretched_duration,30*pi/180,sim.lambda0*1e9,t,prop_output.fields,1e-3/1000,100 );
% Stretched method 2
chirp_sign = 1; % for positive chirp (normal dispersive delay)
omega = ifftshift(2*pi*f,1);
spectrum_amplitude = ifft(prop_output.fields);
func_chirp = calc_chirp;
%[chirp,stretched_field] = func_chirp.General( stretched_duration,omega,spectrum_amplitude,chirp_sign );
prop_output.fields = struct('forward', stretched_field,...
                            'backward',zeros(Nt,1));

DT = 1/intra_burst_rate*1e12; % ps
time_window = zeros(1,num_pulses*2);
time_window(1) = 1/burst_rep_rate*1e12 - DT*num_pulses + (DT-pulse_time_window); % window without pulses
time_window(2) = pulse_time_window; % window with pulses
for i = 2:num_pulses
    time_window(i*2-1) = DT - pulse_time_window; % window without pulses
    time_window(i*2) = pulse_time_window; % window with pulses
end

initial_condition(1,num_pulses) = struct('dt',[],'fields',[]);
for i = 1:num_pulses
    % window without pulses
    initial_condition(i*2-1) = struct('dt',time_window(i*2-1)/Nt,...
                                      'fields',struct('forward', zeros(Nt,1),...
                                                      'backward',zeros(Nt,1)));
    % window with pulses
    initial_condition(i*2) = prop_output;
end

%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

%% Finish the simulation and save the data
save('Yb_CPA_wASE.mat','-v7.3');