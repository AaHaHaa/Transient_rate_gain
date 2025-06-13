% This code solves the 6-um-core Yb-doped gain-managed fiber amplifier 
% (GMNA).

close all; clearvars;

addpath('../../../../../../GMMNLSE/GMMNLSE algorithm/','../../../../../../GMMNLSE/user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_GMMNLSE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1080e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
sim.gpu_yes = false;
sim.save_period = 0.03;
sim.dz = 1e-3;

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain (6um)';
fiber.L0 = 0.3; % m; fiber length
fiber.MFD = 71; % um; mode-field diameter; this is the default settings loaded below

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode'); % load default parameters for "fiber" and "sim"

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% Activating "reuse_data" or "linear_oscillator_model" requires setting other parameters.
% Check the example or "gain_info.m".
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 95; % um
gain_rate_eqn.cladding_diameter = 400; % um
gain_rate_eqn.core_NA = 0.11;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 915; % nm
gain_rate_eqn.absorption_to_get_N_total = 7.5; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 14; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = true; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/10e3; % Assume 15 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^14; % the number of time points
time_window = 500; % ps
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
tfwhm = 0.3; % ps
total_energy = 0.1; % nJ

pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

stretched_duration = 200; % ps
[~,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',stretched_duration,30*pi/180,sim.lambda0*1e9,t,prop_output.fields,1e-3/1000,1 );

prop_output.fields = stretched_field;

%% Run the simulation
num_passes = 6;
num_z = fiber.L0/sim.dz + 1;
kkk = 1;
num_rep = 5;
energy = zeros((num_rep+1)*num_passes,1);
for i = 1:num_passes
    prop_output = GMMNLSE_propagate(fiber, prop_output, sim, gain_rate_eqn);

    saved_data{i} = prop_output.saved_data.signal_fields;
    if i > 1
        for j = 1:num_z
            gain_rate_eqn.saved_data.signal_fields_backward{j} = sqrt(abs(prop_output.saved_data.signal_fields{num_z-j+1}).^2 + abs(gain_rate_eqn.saved_data.signal_fields_backward{j}).^2);
        end
    else
        gain_rate_eqn.saved_data = prop_output.saved_data;
        gain_rate_eqn.saved_data.signal_fields_backward = gain_rate_eqn.saved_data.signal_fields;
    end

    energy(kkk) = trapz(t,abs(prop_output.fields(:,:,end)).^2)/1e3;
    close all;figure;plot(energy(1:kkk));
    kkk = kkk + 1;
end

for ii = 1:num_rep
    for i = 1:num_passes
        v = 1:num_passes; v(i) = [];
        for k = v
            gain_rate_eqn.saved_data.signal_fields_backward = cell(1,1,num_z);
            gain_rate_eqn.saved_data.signal_fields_backward(:) = {zeros(Nt,1)};
            for j = 1:num_z
                gain_rate_eqn.saved_data.signal_fields_backward{j} = sqrt(abs(gain_rate_eqn.saved_data.signal_fields_backward{j}).^2 + abs(saved_data{k}{num_z-j+1}).^2);
            end
        end
    
        if i == 1
            prop_output.fields = stretched_field;
        end
    
        prop_output = GMMNLSE_propagate(fiber, prop_output, sim, gain_rate_eqn);
    
        saved_data{i} = prop_output.saved_data.signal_fields;
    
        energy(kkk) = trapz(t,abs(prop_output.fields(:,:,end)).^2)/1e3;
        close all;figure;plot(energy(1:kkk));
        kkk = kkk + 1;
    end
end

%% Finish the simulation and save the data
%[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,prop_output.fields(:,:,end),'Treacy-t',pi/6,1e-6 );

%func = analyze_sim;
%fig_gain = func.analyze_gain(prop_output.z,[],prop_output.Power.pump,prop_output.population);

%% Plot
figure;
plot(energy/1e3,'linewidth',2,'Color','b');
xlim([1,num_passes*(num_rep+1)]);
xlabel('Iteration');
ylabel('Pulse energy (µJ)');
set(gca,'fontsize',20);
print(gcf,'energy.pdf','-dpdf');

figure;
plot(prop_output.z,squeeze(prop_output.population)*100,'linewidth',2,'Color','b');
xlabel('z (m)');
ylabel('N_1 (%)')
set(gca,'fontsize',20);
print(gcf,'N1.pdf','-dpdf');

clearvars fig_gain;
save('steady_state_regen.mat');