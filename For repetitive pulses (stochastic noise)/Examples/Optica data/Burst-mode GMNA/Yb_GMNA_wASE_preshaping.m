% This code solves the 6-um-core burst-mode GMNA with seed preshaping.
% ASE is included.
% 
% Co-pump power:
% For 4 pulses, use 0.67 W.
%    10 pulses, use 1.1 W.
%    50 pulses, use 4.3 W.
%   100 pulses, use 7.7 W.
%   200 pulses, use  17 W.

close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

num_pulses = 200;

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1080e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
%sim.gpu_yes = false;
sim.save_period = 0.1;
sim.gpuDevice.Index = 1;

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain (6um)';
fiber.L0 = 2.5; % m; fiber length
%fiber.MFD = 6; % um; mode-field diameter; this is the default settings loaded below

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
gain_rate_eqn.copump_power = 17; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 5; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^(nextpow2(2^18/100*num_pulses)); % the number of time points
dt = 0.01;
time_window = Nt*dt; % ps
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

sim.incoherent.dt = dt*2;
incoherent_lambda0 = 1070e-9; % m
sim.incoherent.f0 = 2.99792458e-4/incoherent_lambda0;
incoherent_f = sim.incoherent.f0+(-Nt/2:Nt/2-1)'/(Nt*sim.incoherent.dt); % THz
incoherent_lambda = c./(incoherent_f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
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

tfwhm = 0.5; % ps
total_energy = 0.1*num_pulses; % nJ

pulse_lambda0 = 1025e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
DT = 7; % ps
energy_coeff = linspace(1,3,num_pulses);
prop_output = build_MMgaussian(tfwhm, time_window, total_energy*sum(energy_coeff)/num_pulses, num_pulses, Nt, {'ifft',freq_shift}, sqrt(energy_coeff), (-floor(num_pulses/2):ceil(num_pulses/2)-1)*DT);
prop_output.fields = struct('forward', sum(prop_output.fields,2),...
                            'backward',zeros(Nt,1));

prop_output_woPulse = struct('dt',(1/burst_rep_rate*1e12-time_window)/Nt,...
                             'fields',struct('forward', zeros(Nt,1),...
                                             'backward',zeros(Nt,1)));

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output(2).fields.forward).^2,1),2)*prop_output(2).dt/10^3); % energy in nJ

[grating_spacing,~,dechirped_field] = pulse_compressor_bursts( 'Treacy-t',pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
figure;
plot(t,abs(dechirped_field).^2,'linewidth',2);
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);

%pulse_compressor_animator( 'Treacy-t',grating_spacing*linspace(0.7,1.3,71),pi/6,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6,[],[0,1400],sprintf('dechirping_woASE_%upulses',num_pulses) );

%dechirped_field_single = pulse_compressor_single( 'Treacy-t',0.00255,pi/6,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
%figure;
%plot(t,abs(dechirped_field_single).^2,'linewidth',2);
%xlabel('Time (ps)');
%ylabel('Power (W)');
%set(gca,'fontsize',20);

t_woPulse = (-Nt/2:Nt/2-1)'*prop_output_woPulse.dt; % ps

grating_spacings = linspace(1.2,1.7,201)*1e-3;%grating_spacing*linspace(0.9,1.1,101);
FWHM_mean = zeros(1,length(grating_spacings));
FWHM_std = zeros(1,length(grating_spacings));
I_mean = zeros(1,length(grating_spacings));
I_std = zeros(1,length(grating_spacings));
for i = 1:length(grating_spacings)
    dechirped_field_i = pulse_compressor_single( 'Treacy-t',grating_spacings(i),pi/180*19,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );

    I = abs(dechirped_field_i).^2;
    threshold = max(I)/3;
    [this_max_I,~,FWHM,~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2,'MinPeakDistance',DT*0.9);

    I_mean(i) = mean(this_max_I);
    I_std(i) = std(this_max_I);
    FWHM_mean(i) = mean(FWHM);
    FWHM_std(i) = std(FWHM);
end
[~,min_idx] = min(FWHM_std);
grating_spacing_min = grating_spacings(min_idx);
dechirped_field_min = pulse_compressor_single( 'Treacy-t',grating_spacing_min,pi/6,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
I = abs(dechirped_field_min).^2;
threshold = max(I)/3;
[~,~,FWHM,~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2,'MinPeakDistance',DT*0.9);

save(sprintf('GMNA_wASE_%upulses_preshaping.mat',num_pulses));