% This code solves the 6-um-core burst-mode GMNA.

close all; clearvars;

addpath('../../../UPPE algorithm/','../../../user_helpers/');

num_pulses = 10;

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_UPPE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1080e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
sim.gpu_yes = false;
sim.save_period = 0.1;
sim.gpuDevice.Index = 1;
%sim.adaptive_dz.threshold = 1e-6;
%sim.pulse_centering = false;

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
gain_rate_eqn.copump_power = 1; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

gain_rate_eqn.woPulse = struct('lambda',[1000,1100],'Nt',2^10);

%% Setup general parameters
Nt = 2^13; % the number of time points
dt = 0.01;
pulse_time_window = Nt*dt; % ps; for pulses
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
burst_rep_rate = 1e6; % Hz

tfwhm = 0.5; % ps
seed_energy = 0.1; % nJ

pulse_lambda0 = 1025e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
intra_burst_rate = 60e6; % Hz; intra-burst rate
DT = 1/intra_burst_rate*1e12; % ps

time_window = zeros(1,num_pulses*2);
time_window(1) = 1/burst_rep_rate*1e12 - DT*num_pulses + (DT-pulse_time_window); % window without pulses
time_window(2) = pulse_time_window; % window with pulses
for i = 2:num_pulses
    time_window(i*2-1) = DT - pulse_time_window; % window without pulses
    time_window(i*2) = pulse_time_window; % window with pulses
end

prop_output = build_MMgaussian(tfwhm, pulse_time_window, seed_energy, 1, Nt, {'ifft',freq_shift});
prop_output.fields = struct('forward', prop_output.fields,...
                            'backward',zeros(Nt,1));
for i = 1:num_pulses
    % window without pulses
    initial_condition(i*2-1) = struct('dt',time_window(i*2-1)/gain_rate_eqn.woPulse.Nt,...
                                      'fields',struct('forward', zeros(gain_rate_eqn.woPulse.Nt,1),...
                                                      'backward',zeros(gain_rate_eqn.woPulse.Nt,1)));
    % window with pulses
    initial_condition(i*2) = prop_output;
end


%% Run the simulation
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output(2).fields.forward).^2,1),2)*prop_output(2).dt/10^3); % energy in nJ

% Plot all pulses
% Note the potential plotting gap between time windows
figure;
yyaxis right;
this_t = 0;
for i = 1:num_pulses*2
    if mod(i,2) == 1
        t_tmp = (-gain_rate_eqn.woPulse.Nt/2:gain_rate_eqn.woPulse.Nt/2-1)'*prop_output(i).dt;
        if i == 1
            this_t = t_tmp - t_tmp(1);
        else
            this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        end
        plot(this_t/1e6,zeros(gain_rate_eqn.woPulse.Nt,1),'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    else
        t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
        this_t = (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
        plot(this_t/1e6,abs(prop_output(i).fields.forward(:,:,end)).^2/1e3,'Color',[0.8510,0.3255,0.0980],'LineStyle','-','Marker','None','linewidth',2);
    end
    hold on;
end
hold off;
xlabel('Time (\mus)')
ylabel('Power (kW)');
set(gca,'YColor',[0.8510,0.3255,0.0980]);

yyaxis left;
this_t = 0;
for i = 1:num_pulses*2
    t_tmp = (-size(prop_output(i).population,1)/2:size(prop_output(i).population,1)/2-1)'*prop_output(i).dt;
    if i == 1
        this_t = t_tmp - t_tmp(1);
    else
        this_t =  (this_t(end) + (this_t(end)-this_t(end-1))) + (t_tmp - t_tmp(1));
    end
    plot(this_t/1e6,prop_output(i).population(:,:,end)*100,'Color','b','LineStyle','-','Marker','None','linewidth',2);
    hold on;
end
hold off;
%xlabel('Time (\mus)')
ylabel('N_1 (%)');
set(gca,'YColor','b');
set(gca,'fontsize',20);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
print(gcf,'all_pulses.pdf','-dpdf','-bestfit','-vector');

grating_spacings = linspace(1,1.5,201)*1e-3;%grating_spacing*linspace(0.9,1.1,101);
FWHM_mean = zeros(1,length(grating_spacings));
FWHM_std = zeros(1,length(grating_spacings));
I_mean = zeros(1,length(grating_spacings));
I_std = zeros(1,length(grating_spacings));
for i = 1:length(grating_spacings)
    I_max = zeros(1,num_pulses);
    FWHM = zeros(1,num_pulses);
    for j = 1:num_pulses
        dechirped_field_j = pulse_compressor_single( 'Treacy-t',grating_spacings(i),pi/180*19,sim.lambda0*1e9,t,prop_output(j*2).fields.forward(:,:,end),1e-6 );

        I = abs(dechirped_field_j).^2/1e6;
        threshold = max(I)/1.01;
        [I_max(j),~,FWHM(j),~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
    end
    I_mean(i) = mean(I_max);
    I_std(i) = std(I_max);
    FWHM_mean(i) = mean(FWHM);
    FWHM_std(i) = std(FWHM);
end
[~,min_idx] = min(FWHM_std);
grating_spacing_min = grating_spacings(min_idx);
dechirped_field_min = pulse_compressor_single( 'Treacy-t',grating_spacing_min,pi/6,sim.lambda0*1e9,t,prop_output(2).fields.forward(:,:,end),1e-6 );
I = abs(dechirped_field_min).^2;
threshold = max(I)/1.01;
[~,~,FWHM,~] = findpeaks(I,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);

% Dechirped information
figure;
yyaxis left;
plot(grating_spacings*1e3,I_mean,'linewidth',2','Color','b');
ylabel('avg. peak power (MW)');
set(gca,'YColor','b');
yyaxis right;
plot(grating_spacings*1e3,I_std,'linewidth',2','Color','r');
ylabel('\sigma_{peak power} (MW)');
xlabel('grating spacing (mm)');
set(gca,'YColor','r');
set(gca,'fontsize',20);
xlim([min(grating_spacings),max(grating_spacings)]*1e3);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
print(gcf,'dechirped_info.pdf','-dpdf','-bestfit');

% Spectrum
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain
spectrum = zeros(Nt,num_pulses);
for i = 1:num_pulses
    spectrum(:,i) = abs(fftshift(ifft(prop_output(i*2).fields.forward(:,:,end)),1)).^2*factor_correct_unit.*factor;
end
spectrum_mean = mean(spectrum,2);
spectrum_std = std(spectrum,[],2);
figure;
confplot(lambda,spectrum_mean,spectrum_std,'linewidth',2,'Color','b');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('PSD (nJ/nm)');
xlim([950,1200]);
print(gcf,'spectrum_all.pdf','-dpdf');

figure;
confplot(lambda,spectrum_mean,spectrum_std,'linewidth',2,'Color','b');
set(gca,'fontsize',30);
xlabel('Wavelength (nm)');
xlim([1060,1100]);
ylim([0.335,0.355]);
print(gcf,'spectrum_all (closeview).pdf','-dpdf');

% Temporal profiles of the first pulse
figure;
yyaxis right;
plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2,'Color',[0.8510,0.3255,0.0980]);
xlabel('Time (ps)');
ylabel('Power (kW)');
set(gca,'YColor',[0.8510,0.3255,0.0980]);
set(gca,'fontsize',20);
xlim([-5,5]);
yyaxis left;
plot(t,prop_output(2).population(:,:,end)*100,'linewidth',2,'Color','b');
set(gca,'YColor','b');
ylabel('N_1 (%)');
YLim = get(gca,'YLim'); ylim([min(YLim)*0.9999,max(YLim)*1.0001]);
print(gcf,'temporal profile1.pdf','-dpdf');

% Spectrum of the first pulse
figure;
plot(lambda,spectrum(:,1),'linewidth',2,'Color','b');
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('PSD (nJ/nm)');
xlim([950,1200]);
print(gcf,'spectrum1.pdf','-dpdf');

save(sprintf('GMNA_woASE_%upulses.mat',num_pulses));