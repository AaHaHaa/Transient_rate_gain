% This code runs a fiber-based regen with Coherent's XLMA-YTF-100/400/480
% fiber.
%
% Regen has a 60-MHz cavity repetition rate, and 10-kHz operation rate.
%
% Reference:
% Haig et al., "Single-mode regenerative amplification in multimode fiber,"
% Optica 10(11), 1417-1420 (2024).

close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

num_pulses = 6;

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1030e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
sim.gpu_yes = false;
sim.save_period = 0.1;
%sim.gpuDevice.Index = 1;
%sim.pulse_centering = false;

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain (6um)';
fiber.L0 = 0.3; % m; fiber length
fiber.MFD = 71; % um; mode-field diameter; this is the default settings loaded below

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
gain_rate_eqn.core_diameter = 95; % um
gain_rate_eqn.cladding_diameter = 400; % um
gain_rate_eqn.core_NA = 0.11;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 915; % nm
gain_rate_eqn.absorption_to_get_N_total = 7.5; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 14; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^14; % the number of time points
dt = 0.05;
pulse_time_window = Nt*dt; % ps; for pulses
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
burst_rep_rate = 10e3; % Hz

tfwhm = 0.3; % ps
seed_energy = 0.1; % nJ

pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
intra_burst_rate = 70e6; % Hz; intra-burst rate
DT = 1/intra_burst_rate*1e12; % ps

time_window = zeros(1,num_pulses*2);
time_window(1) = 1/burst_rep_rate*1e12 - DT*num_pulses + (DT-pulse_time_window); % window without pulses
time_window(2) = pulse_time_window; % window with pulses
for i = 2:num_pulses
    time_window(i*2-1) = DT - pulse_time_window; % window without pulses
    time_window(i*2) = pulse_time_window; % window with pulses
end

prop_output = build_MMgaussian(tfwhm, pulse_time_window, seed_energy, 1, Nt, {'ifft',freq_shift});
stretched_duration = 200; % ps
[~,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',stretched_duration,30*pi/180,sim.lambda0*1e9,t,prop_output.fields,1e-3/1000,1 );
prop_output.fields = struct('forward', stretched_field,...
                            'backward',zeros(Nt,1));
for i = 1:num_pulses
    % window without pulses
    initial_condition(i*2-1) = struct('dt',time_window(i*2-1)/Nt,...
                                      'fields',struct('forward', zeros(Nt,1),...
                                                      'backward',zeros(Nt,1)));
    % window with pulses
    initial_condition(i*2) = prop_output;
end


%% Run the simulation
prop_output = initial_condition;
last_rt = 20;
last_energy = zeros(last_rt,1);
for i = 1:last_rt
    prop_output = Transient_gain_UPPE_propagate(fiber, prop_output, sim, gain_rate_eqn);

    last_energy(i) = sum(abs(prop_output(end).fields.forward(:,:,end)).^2)*prop_output(end).dt/1e3;
    disp(last_energy(i));

    if i < last_rt
        for window_i = num_pulses*2:-2:4
            prop_output(window_i).fields.forward = prop_output(window_i-2).fields.forward(:,:,end);
        end
        prop_output(2).fields.forward = initial_condition(2).fields.forward;
    
        for window_i = 1:2:num_pulses*2
            prop_output(window_i).fields.forward = zeros(size(initial_condition(window_i).fields.forward));
        end
    end
end

%% Finish the simulation and save the data
% Energy
figure;
plot(last_energy/1e3,'linewidth',2,'Color','b');
xlabel('Iteration');
ylabel('Last pulse''s energy (µJ)');
set(gca,'fontsize',20);
print(gcf,'last_energy.pdf','-dpdf');

energy = zeros(num_pulses,1);
for i = 1:num_pulses
    energy(i) = sum(abs(prop_output(i*2).fields.forward(:,:,end)).^2)*prop_output(i*2).dt/1e3;
end
figure;
plot(energy,'linewidth',2,'Color','b');
xlabel('Pulse');
ylabel('Energy (nJ)');
set(gca,'fontsize',20);
print(gcf,'energy.pdf','-dpdf');

%% Plot all pulses
% Note the potential plotting gap between time windows
figure;
yyaxis right;
this_t = 0;
for i = 1:num_pulses*2
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
for i = 1:num_pulses*2
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
%xlabel('Time (ps)')
ylabel('N_1 (%)');
set(gca,'YColor','b');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
%print(gcf,'pulse_t.pdf','-dpdf','-bestfit');
exportgraphics(gcf, 'pulse_t.pdf', 'ContentType', 'vector');

xlim([time_window(1),this_t(end)]/1e6);
%print(gcf,'pulse_t(zoomIn).pdf','-dpdf','-bestfit');
exportgraphics(gcf, 'pulse_t(zoomIn).pdf', 'ContentType', 'vector');

%% Plot all pulses (log scale)
% Note the potential plotting gap between time windows
figure;
yyaxis right;
this_t = 0;
for i = 1:num_pulses*2
    if mod(i,2) == 1
        t_tmp = (-Nt/2:Nt/2-1)'*prop_output(i).dt;
        if i == 1
            this_t = t_tmp - t_tmp(1);
        else
            this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        end
        semilogy(this_t/1e6,zeros(Nt,1),'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    else
        t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
        last_t_end = this_t(end);
        this_t = (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        semilogy([last_t_end;this_t]/1e6,[0;abs(prop_output(i).fields.forward(:,:,end)).^2/1e3],'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    end
    hold on;
end
hold off;
xlabel('Time (µs)')
ylabel('Power (kW)');
set(gca,'YColor',[0.8510,0.3255,0.0980]);
ylim([0.001,500]);
set(gca,'YTick',[1e-3,1e-1,1e1,1e3]);

yyaxis left;
this_t = 0;
for i = 1:num_pulses*2
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
%xlabel('Time (ps)')
ylabel('N_1 (%)');
set(gca,'YColor','b');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
xlim([time_window(1),this_t(end)]/1e6);
%print(gcf,'pulse_t_log(zoomIn).pdf','-dpdf','-bestfit');
exportgraphics(gcf, 'pulse_t_log(zoomIn).pdf', 'ContentType', 'vector');

%% Spectrum
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain
spectrum = zeros(Nt,num_pulses);
for i = 1:num_pulses
    spectrum(:,i) = abs(fftshift(ifft(prop_output(i*2).fields.forward(:,:,end)),1)).^2*factor_correct_unit.*factor;
end
figure;
plot(lambda,spectrum(:,:,end),'linewidth',2,'Color','b');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('PSD (nJ/nm)');

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power,fig] = analyze_field( t,f,prop_output(end).fields.forward(:,:,end),'Treacy-t',pi/6,1e-6,true );
figure(fig(2));
set(gca,'fontsize',30);
title('');
legend('TL','D');

save('regen.mat');