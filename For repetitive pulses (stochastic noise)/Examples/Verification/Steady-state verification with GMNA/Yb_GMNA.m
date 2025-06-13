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
sim.save_period = 0.1;

% -------------------------------------------------------------------------

% Gain fiber
sim.progress_bar_name = 'Gain';
fiber.L0 = 3; % m; fiber length
fiber.MFD = 5.95; % um; mode-field diameter

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
gain_rate_eqn.copump_power = 1.6; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^14; % the number of time points
time_window = 100; % ps
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
rep_rate = 30e6; % Hz

tfwhm = 0.5; % ps
total_energy = 0.1; % nJ

pulse_lambda0 = 1025e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

prop_output.fields = struct('forward', prop_output.fields,...
                            'backward',zeros(Nt,1));

prop_output_woPulse = struct('dt',(1/rep_rate*1e12-time_window)/Nt,...
                             'fields',struct('forward', zeros(Nt,1),...
                                             'backward',zeros(Nt,1)));

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%% Run the simulation
prop_output = Periodic_transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

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

t_woPulse = (-Nt/2:Nt/2-1)'*prop_output_woPulse.dt; % ps

% Plot all pulses
% Note the potential plotting gap between time windows, so plotting below
% requires care to fill it.
figure;
yyaxis right;
this_t = 0;
for i = 1:2
    if mod(i,2) == 1
        t_tmp = (-Nt/2:Nt/2-1)'*prop_output(i).dt;
        if i == 1
            this_t = t_tmp - t_tmp(1);
        else
            this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        end
        plot(this_t/1e6,zeros(Nt,1),'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    else
        t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
        last_t_end = this_t(end);
        this_t = (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        plot([last_t_end;this_t]/1e6,[0;abs(prop_output(i).fields.forward(:,:,end)).^2/1e3],'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    end
    hold on;
end
hold off;
xlabel('Time (µs)')
ylabel('Power (kW)');
set(gca,'YColor',[0.8510,0.3255,0.0980]);

yyaxis left;
this_t = 0;
for i = 1:2
    t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
    last_t_end = this_t(end);
    if i == 1
        this_t = t_tmp - t_tmp(1);
    else
        this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
    end
    if mod(i,2) == 1
        plot(this_t/1e6,prop_output(i).population(:,:,end)*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
    else
        plot([last_t_end;this_t]/1e6,[prop_output(i-1).population(end,:,end);prop_output(i).population(:,:,end)]*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
    end
    hold on;
end
hold off;
%xlabel('Time (\mus)')
ylabel('N_1 (%)');
set(gca,'YColor','b');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
xlim([0,1/rep_rate*1e6]);
print(gcf,'all_pulses.pdf','-dpdf','-bestfit','-vector');

%% Save results
save('transient_GMNA.mat','f','t','prop_output','sim','fiber','gain_rate_eqn');