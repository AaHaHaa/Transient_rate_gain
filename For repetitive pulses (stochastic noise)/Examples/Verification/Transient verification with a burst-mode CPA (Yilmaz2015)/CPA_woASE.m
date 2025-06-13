close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1034e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
%sim.gpu_yes = false;
sim.save_period = 0.05;

% -------------------------------------------------------------------------

% Gain fiber
sim.progress_bar_name = 'Gain';
fiber.L0 = 2; % m; fiber length
fiber.MFD = 2500; % um; mode-field diameter

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
gain_rate_eqn.absorption_wavelength_to_get_N_total = 976; % nm
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

gain_rate_eqn.woPulse = struct('lambda',[1000,1100],'Nt',2^10);

% Each type
gain_rate_eqn1 = gain_rate_eqn;
gain_rate_eqn1.core_diameter = 6; % um
gain_rate_eqn1.cladding_diameter = 6; % um
gain_rate_eqn1.core_NA = 0.12;
gain_rate_eqn1.absorption_to_get_N_total = 500; % dB/m
gain_rate_eqn1.copump_power = 0.017; % W

gain_rate_eqn2 = gain_rate_eqn;
gain_rate_eqn2.core_diameter = 6; % um
gain_rate_eqn2.cladding_diameter = 6; % um
gain_rate_eqn2.core_NA = 0.12;
gain_rate_eqn2.absorption_to_get_N_total = 1200; % dB/m
gain_rate_eqn2.copump_power = 0.04; % W

gain_rate_eqn3 = gain_rate_eqn;
gain_rate_eqn3.core_diameter = 20; % um
gain_rate_eqn3.cladding_diameter = 125; % um
gain_rate_eqn3.core_NA = 0.08;
gain_rate_eqn3.absorption_to_get_N_total = 29; % dB/m
gain_rate_eqn3.copump_power = 2; % W

gain_rate_eqn4 = gain_rate_eqn;
gain_rate_eqn4.core_diameter = 25; % um
gain_rate_eqn4.cladding_diameter = 250; % um
gain_rate_eqn4.core_NA = 0.06;
gain_rate_eqn4.absorption_to_get_N_total = 6.5; % dB/m
gain_rate_eqn4.copump_power = 30; % W

gain_rate_eqn5 = gain_rate_eqn;
gain_rate_eqn5.core_diameter = 25; % um
gain_rate_eqn5.cladding_diameter = 250; % um
gain_rate_eqn5.core_NA = 0.06;
gain_rate_eqn5.absorption_to_get_N_total = 6.5; % dB/m
gain_rate_eqn5.copump_power = 30; % W

gain_rate_eqn6 = gain_rate_eqn;
gain_rate_eqn6.core_diameter = 25; % um
gain_rate_eqn6.cladding_diameter = 250; % um
gain_rate_eqn6.core_NA = 0.06;
gain_rate_eqn6.absorption_to_get_N_total = 10.8; % dB/m
gain_rate_eqn6.copump_power = 30; % W

%% Setup general parameters
Nt = 2^18; % the number of time points
time_window = 300e3; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn1 = gain_info( fiber,sim,gain_rate_eqn1,ifftshift(lambda,1) );
gain_rate_eqn2 = gain_info( fiber,sim,gain_rate_eqn2,ifftshift(lambda,1) );
gain_rate_eqn3 = gain_info( fiber,sim,gain_rate_eqn3,ifftshift(lambda,1) );
gain_rate_eqn4 = gain_info( fiber,sim,gain_rate_eqn4,ifftshift(lambda,1) );
gain_rate_eqn5 = gain_info( fiber,sim,gain_rate_eqn5,ifftshift(lambda,1) );
gain_rate_eqn6 = gain_info( fiber,sim,gain_rate_eqn6,ifftshift(lambda,1) );

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

%% Fiber length
fiber1 = fiber;
fiber1.L0 = 0.4;

fiber2 = fiber;
fiber2.L0 = 0.15;

fiber3 = fiber;
fiber3.L0 = 1.4;

fiber4 = fiber;
fiber4.L0 = 2.8;

fiber5 = fiber;
fiber5.L0 = 1;

fiber6 = fiber;
fiber6.L0 = 2.7;

%% Setup initial conditions
burst_rep_rate = 1e6; % Hz

num_pulses = 10;

tfwhm = 10; % ps
total_energy = 0.13*num_pulses; % nJ

pulse_lambda0 = 1034e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
DT = 10e3; % ps
energy_coeff = linspace(1,1,num_pulses);
prop_output = build_MMgaussian(tfwhm, time_window, total_energy*sum(energy_coeff)/num_pulses, num_pulses, Nt, {'ifft',freq_shift}, sqrt(energy_coeff), (-floor(num_pulses/2):ceil(num_pulses/2)-1)*DT);
prop_output.fields = struct('forward', sum(prop_output.fields,2),...
                            'backward',zeros(Nt,1));

prop_output_woPulse = struct('dt',(1/burst_rep_rate*1e12-time_window)/gain_rate_eqn.woPulse.Nt,...
                             'fields',struct('forward', zeros(gain_rate_eqn.woPulse.Nt,1),...
                                             'backward',zeros(gain_rate_eqn.woPulse.Nt,1)));

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber1, initial_condition, sim, gain_rate_eqn1);
prop_output(1).fields = struct('forward', zeros(gain_rate_eqn.woPulse.Nt,1),...
                               'backward',zeros(gain_rate_eqn.woPulse.Nt,1));
prop_output(2).fields = struct('forward',prop_output(2).fields.forward(:,:,end),...
                               'backward',zeros(Nt,1));
prop_output = Transient_gain_UPPE_propagate(fiber2, prop_output,       sim, gain_rate_eqn2);
prop_output(1).fields = struct('forward', zeros(gain_rate_eqn.woPulse.Nt,1),...
                               'backward',zeros(gain_rate_eqn.woPulse.Nt,1));
prop_output(2).fields = struct('forward',prop_output(2).fields.forward(:,:,end),...
                               'backward',zeros(Nt,1));
prop_output = Transient_gain_UPPE_propagate(fiber3, prop_output,       sim, gain_rate_eqn3);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output(2).fields.forward).^2,1),2)*prop_output(2).dt/10^3); % energy in nJ

t1 = (-gain_rate_eqn.woPulse.Nt/2:gain_rate_eqn.woPulse.Nt/2-1)'*prop_output_woPulse.dt; % ps

figure;
h = plot(t1/1e6,prop_output(1).population(:,:,end)*100,'linewidth',2);
set(h,'Color','b');
set(gca,'fontsize',20);
xlabel('Time (\mus)'); ylabel('N_2 (%)');
%print(gcf,'N2_1_woASE.jpg','-djpeg');

figure;
yyaxis left;
plot(t,prop_output(2).population(:,:,end)*100,'linewidth',2);
ylabel('N_2 (%)');
yyaxis right;
plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
%print(gcf,'Aout_woASE.jpg','-djpeg');

figure;
plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2,'Color','b');
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);
%print(gcf,'Ain_woASE.jpg','-djpeg');

save(sprintf('Yilmaz_%ukHz.mat',burst_rep_rate/1e3));