close all; clearvars;

addpath('../../UPPE algorithm/','../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1030e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
sim.gpu_yes = false;
sim.pulse_centering = false;
fiber.n2 = 0;

% -------------------------------------------------------------------------

% Gain fiber
sim.progress_bar_name = 'Gain (6um)';
fiber.L0 = 10; % m; fiber length
sim.dz = 0.01;
save_num = 50; % the number of saved data
sim.save_period = fiber.L0/save_num;

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
gain_rate_eqn.copump_power = 0.6; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 20; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information (final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^11; % the number of time points
time_window = 5; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( sim,gain_rate_eqn,{ifftshift(lambda,1),ifftshift(lambda,1)} );

sim.incoherent_dt = dt;

%% calculate fiber betas from silica refractive index
% This is important to correctly simulate the broadband situations.
% Taylor-series coefficients is only good in narrowband situations.

% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9)*0;

%% Setup initial conditions
burst_rep_rate = 10e3; % Hz

num_pulses = 1;

tfwhm = 0.1; % ps
total_energy = 1e-100*num_pulses; % nJ

pulse_lambda0 = 1025e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
DT = 3; % ps
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, num_pulses, Nt, {'ifft',freq_shift}, ones(1,num_pulses), (-floor(num_pulses/2):ceil(num_pulses/2)-1)*DT);
prop_output.fields = struct('forward', sum(prop_output.fields,2),...
                            'backward',zeros(Nt,1));

prop_output_woPulse = struct('dt',(1/burst_rep_rate*1e12-time_window)/Nt,...
                             'fields',struct('forward', zeros(Nt,1),...
                                             'backward',zeros(Nt,1)));

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

%% Plot results
ASE = struct('forward', squeeze(mean(abs(prop_output(1).fields.forward).^2,1))*1e3,... % mW
             'backward',squeeze(mean(abs(prop_output(1).fields.backward).^2,1))*1e3);  % mW

figure;
h = plot(prop_output(1).z,[ASE.forward,ASE.backward],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r');
set(gca,'fontsize',20);
xlim([prop_output(1).z(1),prop_output(1).z(end)]);
xlabel('Propagation length (m)'); ylabel('ASE power (mW)');
%print(gcf,'transient ASE.jpg','-djpeg');

t1 = (-Nt/2:Nt/2-1)'*prop_output_woPulse.dt; % ps

figure;
ASE_forward = abs(prop_output(1).fields.forward(:,:,end)).^2*1e3;
avg_ASE_forward = smooth(ASE_forward,30);
h = plot(t1/1e3,[ASE_forward,avg_ASE_forward,ASE.forward(end)*ones(Nt,1)],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','k');
set(gca,'fontsize',20);
xlabel('Time (\mus)'); ylabel('Power (mW)');
legend('|A(t)|^2','smoothed');
%print(gcf,'ASE_forward(t).jpg','-djpeg');

figure;
ASE_backward = abs(prop_output(1).fields.backward(:,:,1)).^2*1e3;
avg_ASE_backward = smooth(ASE_backward,30);
h = plot(t1/1e3,[ASE_backward,avg_ASE_backward,ASE.backward(1)*ones(Nt,1)],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','k');
set(gca,'fontsize',20);
xlabel('Time (\mus)'); ylabel('Power (mW)');
legend('|A(t)|^2','smoothed');
%print(gcf,'ASE_backward(t).jpg','-djpeg');

%%
figure;
yyaxis left;
plot(t1,prop_output(1).population(:,:,end)*100,'linewidth',2);
ylabel('N_1 (%)');
yyaxis right;
plot(t1,abs(prop_output(1).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);

figure;
yyaxis left;
plot(t,prop_output(2).population(:,:,end)*100,'linewidth',2);
ylabel('N_1 (%)');
yyaxis right;
plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
xlabel('Time (ps)'); ylabel('Power (kW)');
set(gca,'fontsize',20);

%% Save results
save('transient_ASE_evolutions.mat','f','t','prop_output','sim','fiber','gain_rate_eqn');