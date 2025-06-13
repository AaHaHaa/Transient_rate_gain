close all; clearvars;

addpath('../../../user_helpers/');

%% Load general information
load('oscillator_general_info.mat');

f1_real = linspace(c/gain_rate_eqn.woPulse.lambda(2),c/gain_rate_eqn.woPulse.lambda(1),gain_rate_eqn.woPulse.Nt)'; % THz
lambda1 = c./f1_real; % nm
factor_correct_unit1 = (gain_rate_eqn.woPulse.Nt*dt_woPulse)^2/1e3; % nJ/THz
factor1 = c./lambda1.^2; % change the spectrum from frequency domain into wavelength domain
f1 = (-gain_rate_eqn.woPulse.Nt/2:gain_rate_eqn.woPulse.Nt/2-1)'/(gain_rate_eqn.woPulse.Nt*dt_woPulse); % THz
correction_factor1 = (max(f1)-min(f1))/(max(f1_real)-min(f1_real));

factor_correct_unit2 = (Nt*dt)^2/1e3; % nJ/THz
factor2 = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

%% Evolutions (in the supplement)
total_rt = 5000;
rt = (1:total_rt)';

energy_all = zeros(total_rt,1);
peak_power_all = zeros(total_rt,1);
population1_all = zeros(total_rt,1);
population2_all = zeros(total_rt,1);
population3_all = zeros(total_rt,1);
population4_all = zeros(total_rt,1);
dlambda1_all = zeros(total_rt,1);

name_idx = 1000:1000:total_rt;
for i = 1:length(name_idx)
    load(sprintf('oscillator_%u.mat',name_idx(i)));

    energy_all(1000*(i-1)+1:1000*i) = sum(output_energy,2);
    peak_power_all(1000*(i-1)+1:1000*i) = max([squeeze(max(abs(output_field{1}).^2,[],1)),squeeze(max(abs(output_field{2}).^2,[],1))],[],2);
    population1_all(1000*(i-1)+1:1000*i) = population1;
    population2_all(1000*(i-1)+1:1000*i) = population2;
    population3_all(1000*(i-1)+1:1000*i) = population3;
    population4_all(1000*(i-1)+1:1000*i) = population4;

    spectrum1 = abs(fftshift(ifft(output_field{1}),1)).^2*factor_correct_unit1; % nJ/THz
    spectrum1 = spectrum1.*factor1*correction_factor1; % nJ/nm
    spectrum1 = spectrum1*1e-9/(gain_rate_eqn.woPulse.Nt*dt_woPulse*1e-12); % W/nm
    spectrum1 = interp1(lambda1,spectrum1,lambda,'linear',0);
    spectrum2 = abs(fftshift(ifft(output_field{2}),1)).^2*factor_correct_unit2; % nJ/THz
    spectrum2 = spectrum2.*factor2; % nJ/nm
    spectrum2 = spectrum2*1e-9/(Nt*dt*1e-12); % W/nm
    dlambda1_all(1000*(i-1)+1:1000*i) = calc_RMS(lambda,spectrum1+spectrum2);
end
%{
figure;
plot(rt/1e3,dlambda1_all,'linewidth',2,'Color','b');
xlim([2,3]);
xlabel('Roundtrip (x10^3)');
ylabel('RMS bandwidth (nm)');
ax = gca;
set(ax,'fontsize',20);
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]); % make the figure "wider"
pos = ax.Position; set(ax, 'Position', [0.143,pos(2),0.77,pos(4)]); % make the axes more narrower
exportgraphics(gcf, 'df_evolution_zoomIn.pdf', 'ContentType', 'vector'); % print command doesn't work for figures with a far-away ylabel

%% Plot single spectrum
load('oscillator_3000.mat');

all_idx = [107,120,130,148,164,400,1000];

for i = 1:length(all_idx)
    idx = all_idx(i);

    spectrum1 = abs(fftshift(ifft(output_field{1}(:,:,idx)),1)).^2*factor_correct_unit1; % nJ/THz
    spectrum1 = spectrum1.*factor1*correction_factor1; % nJ/nm
    spectrum1 = spectrum1*1e-9/(gain_rate_eqn.woPulse.Nt*dt_woPulse*1e-12); % W/nm
    spectrum1 = interp1(lambda1,spectrum1,lambda,'linear',0);
    spectrum2 = abs(fftshift(ifft(output_field{2}(:,:,idx)),1)).^2*factor_correct_unit2; % nJ/THz
    spectrum2 = spectrum2.*factor2; % nJ/nm
    spectrum2 = spectrum2*1e-9/(Nt*dt*1e-12); % W/nm
    spectrum = spectrum1 + spectrum2;
    spectrum = spectrum./max(spectrum);

    figure;
    plot(lambda,spectrum,'linewidth',6,'Color','b');
    xlabel('Wavelength (nm)');
    ylabel('PSD (norm.)');
    set(gca,'fontsize',40);
    xlim([1000,1060]);
    print(gcf,sprintf('f_modelocking%uRT.pdf',idx+2000),'-dpdf');

    time_profile = abs(output_field{1}(:,:,idx)).^2;
    time_profile = time_profile/max(time_profile);
    figure;
    plot(t1,time_profile,'linewidth',6,'Color',[0.851,0.3255,0.098]);
    xlabel('Time (ps)');
    ylabel('Power (norm.)');
    set(gca,'fontsize',40);

    time_profile = abs(output_field{2}(:,:,idx)).^2;
    time_profile = time_profile/max(time_profile);
    figure;
    plot(t,time_profile,'linewidth',6,'Color',[0.851,0.3255,0.098]);
    xlabel('Time (ps)');
    ylabel('Power (norm.)');
    set(gca,'fontsize',40);
    print(gcf,sprintf('t_modelocking%uRT.pdf',idx+2000),'-dpdf');
end
%}
%% Q-switching
figure;
plot(rt/1e3,dlambda1_all,'linewidth',2,'Color','b');
xlim([1.6,2]);
xlabel('Roundtrip (x10^3)');
ylabel('RMS bandwidth (nm)');
ax = gca;
set(ax,'fontsize',20);
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]); % make the figure "wider"
pos = ax.Position; set(ax, 'Position', [0.143,pos(2),0.77,pos(4)]); % make the axes more narrower
exportgraphics(gcf, 'df_Qevolution_zoomIn.pdf', 'ContentType', 'vector'); % print command doesn't work for figures with a far-away ylabel

%% Plot single spectrum

load('oscillator_2000.mat');

all_idx = [600,650,712,727,742,770,800];

for i = 1:length(all_idx)
    idx = all_idx(i);

    spectrum1 = abs(fftshift(ifft(output_field{1}(:,:,idx)),1)).^2*factor_correct_unit1; % nJ/THz
    spectrum1 = spectrum1.*factor1*correction_factor1; % nJ/nm
    spectrum1 = spectrum1*1e-9/(gain_rate_eqn.woPulse.Nt*dt_woPulse*1e-12); % W/nm
    spectrum1 = interp1(lambda1,spectrum1,lambda,'linear',0);
    spectrum2 = abs(fftshift(ifft(output_field{2}(:,:,idx)),1)).^2*factor_correct_unit2; % nJ/THz
    spectrum2 = spectrum2.*factor2; % nJ/nm
    spectrum2 = spectrum2*1e-9/(Nt*dt*1e-12); % W/nm
    spectrum = spectrum1 + spectrum2;
    spectrum = spectrum./max(spectrum);

    figure;
    plot(lambda,spectrum,'linewidth',6,'Color','b');
    xlabel('Wavelength (nm)');
    ylabel('PSD (norm.)');
    set(gca,'fontsize',40);
    xlim([1000,1060]);
    print(gcf,sprintf('f_Q%uRT.pdf',idx+1000),'-dpdf');

    time_profile = abs(output_field{1}(:,:,idx)).^2;
    time_profile = time_profile/max(time_profile);
    figure;
    plot(t1/1e3,time_profile,'linewidth',6,'Color',[0.851,0.3255,0.098]);
    xlabel('Time (ns)');
    ylabel('Power (norm.)');
    set(gca,'fontsize',40);
    print(gcf,sprintf('t_Q%uRT.pdf',idx+1000),'-dpdf');
end