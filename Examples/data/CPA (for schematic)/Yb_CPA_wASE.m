% This code solves the 6-um-core Yb-doped gain-managed fiber amplifier 
% (GMNA).

close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1030e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
%sim.gpu_yes = false;
sim.save_period = 0;
sim.dz = 1e-2;

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
[fiber,sim] = load_default_UPPE_propagate(fiber,sim,'single_mode'); % load default parameters for "fiber" and "sim"

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
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

gain_rate_eqn.woPulse = struct('lambda',[1000,1100],'Nt',2^10);

%% Setup general parameters
Nt = 2^16; % the number of time points
time_window = 4e3; % ps
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
rep_rate = 100e3; % Hz

tfwhm = 0.4; % ps
total_energy = 0.1; % nJ

pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

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

prop_output_woPulse = struct('dt',(1/rep_rate*1e12-time_window)/gain_rate_eqn.woPulse.Nt,...
                             'fields',struct('forward', zeros(gain_rate_eqn.woPulse.Nt,1),...
                                             'backward',zeros(gain_rate_eqn.woPulse.Nt,1)));

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output(2).fields.forward).^2,1),2)*prop_output(2).dt/10^3); % energy in nJ

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power,fig] = analyze_field( t,f,prop_output(2).fields.forward(:,:,end),'Treacy-t',pi/6,1e-6,true );

[grating_spacing,~,dechirped_field] = pulse_compressor_bursts( 'Treacy-t',pi/6,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
figure;
plot(t,abs(dechirped_field).^2,'linewidth',2);
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);

t1 = (-gain_rate_eqn.woPulse.Nt/2:gain_rate_eqn.woPulse.Nt/2-1)'*prop_output_woPulse.dt; % ps

figure;
h = plot(t1/1e6,prop_output(1).population(:,:,end)*100,'linewidth',2);
set(h,'Color','b');
set(gca,'fontsize',20);
xlabel('Time (ms)'); ylabel('N_1 (%)');
print(gcf,'N2_1_wASE.jpg','-djpeg');

figure;
yyaxis left;
plot(t1,prop_output(1).population(:,:,end)*100,'linewidth',2,'Color','b');
ylabel('N_1 (%)');
yyaxis right;
plot(t1,abs(prop_output(1).fields.forward(:,:,end)).^2/1e3 + rand(gain_rate_eqn.woPulse.Nt,1)*1e4,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
%print(gcf,'Aout_1_wASE.jpg','-djpeg');
print(gcf,'Aout_1_wASE.pdf','-dpdf','-bestfit');

figure;
yyaxis left;
plot(t,prop_output(2).population(:,:,end)*100,'linewidth',2,'Color','b');
ylabel('N_1 (%)');
yyaxis right;
plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
%print(gcf,'Aout_2_wASE.jpg','-djpeg');
print(gcf,'Aout_2_wASE.pdf','-dpdf');

figure;
plot(t,abs(prop_output(2).fields.forward(:,:,1)).^2/1e3,'linewidth',2,'Color','b');
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
print(gcf,'Ain_wASE.jpg','-djpeg');

close all;
save('Yb_CPA_wASE.mat');