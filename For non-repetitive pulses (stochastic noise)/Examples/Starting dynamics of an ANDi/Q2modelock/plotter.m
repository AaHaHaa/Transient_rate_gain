close all; clearvars;

addpath('../../../user_helpers/');

%% Load general info
load('oscillator_general_info.mat');

f1_real = linspace(c/gain_rate_eqn.woPulse.lambda(2),c/gain_rate_eqn.woPulse.lambda(1),gain_rate_eqn.woPulse.Nt)'; % THz
lambda1 = c./f1_real; % nm
factor_correct_unit1 = (gain_rate_eqn.woPulse.Nt*dt_woPulse)^2/1e3; % nJ/THz
factor1 = c./lambda1.^2; % change the spectrum from frequency domain into wavelength domain
f1 = (-gain_rate_eqn.woPulse.Nt/2:gain_rate_eqn.woPulse.Nt/2-1)'/(gain_rate_eqn.woPulse.Nt*dt_woPulse); % THz
correction_factor1 = (max(f1)-min(f1))/(max(f1_real)-min(f1_real));

factor_correct_unit2 = (Nt*dt)^2/1e3; % nJ/THz
factor2 = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

%% Evolutions (in the manuscript)
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

fig = figure;
ax1 = axes();
ax2 = copyobj(ax1,fig); % create two identical overlaying axes
ax3 = copyobj(ax1,fig); % create two identical overlaying axes
linkaxes([ax1,ax2,ax3],'x'); % link the x axes
% Energy and peak power
plot(ax1,rt/1e3,energy_all/1e3,'Color','b','linewidth',4); ylabel(ax1,'Energy (\muJ)'); ylim(ax1,[0,3]);
plot(ax2,rt/1e3,peak_power_all/1e3,'Color','r','linewidth',4); ylabel(ax2,'Peak power (kW)'); ylim(ax2,[0,90]);
% Populations
h = plot(ax3,rt/1e3,[population1_all,population2_all,population3_all,population4_all]*100,'linewidth',4); ylabel(ax3,'N_1 (%)'); ylim(ax3,[0,55]); set(ax3, 'YAxisLocation', 'right'); % set y-axes to right side
set(h(1),'Color',0*ones(1,3)); set(h(2),'Color',0.4*ones(1,3)); set(h(3),'Color',0.6*ones(1,3)); set(h(4),'Color',0.7*ones(1,3));
xlabel(ax1,'Roundtrip (x10^3)'); xlabel(ax2,''); xlabel(ax3,'');
set([ax1,ax2,ax3],'Color', 'None','fontsize',20);
set(ax1,'YColor','b'); set(ax2,'YColor','r'); set(ax3,'YColor','k');
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]); % make the figure "wider"
pos = ax1.Position; set([ax1,ax2,ax3], 'Position', [0.143,pos(2),0.77,pos(4)]); % make the axes more narrower
% Extend the ax1 y-axes ticklabels rightward a bit by adding space
ax1.YTickLabel = strcat(ax1.YTickLabel,{'            '});
%print(gcf,'evolution.pdf','-dpdf','-bestfit');
exportgraphics(gcf, 'evolution.pdf', 'ContentType', 'vector'); % print command doesn't work for figures with a far-away ylabel

%% Evolutions (spectrum)
figure;
plot(rt/1e3,dlambda1_all,'linewidth',2,'Color','b');
xlabel('Roundtrip (x10^3)');
ylabel('RMS bandwidth (nm)');
ax = gca;
set(ax,'fontsize',20);
pos = get(gcf,'Position'); set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]); % make the figure "wider"
pos = ax.Position; set(ax, 'Position', [0.143,pos(2),0.77,pos(4)]); % make the axes more narrower
exportgraphics(gcf, 'df_evolution.pdf', 'ContentType', 'vector'); % print command doesn't work for figures with a far-away ylabel

load(sprintf('oscillator_%u.mat',1000));
% Q-switching
spectrum1_588 = abs(fftshift(ifft(output_field{1}(:,:,588)),1)).^2*factor_correct_unit1; % nJ/THz
spectrum1_588 = spectrum1_588.*factor1*correction_factor1; % nJ/nm
spectrum1_588 = spectrum1_588*1e-9/(gain_rate_eqn.woPulse.Nt*dt_woPulse*1e-12); % W/nm
spectrum1_588 = interp1(lambda1,spectrum1_588,lambda,'linear',0);
spectrum2_588 = abs(fftshift(ifft(output_field{2}(:,:,588)),1)).^2*factor_correct_unit2; % nJ/THz
spectrum2_588 = spectrum2_588.*factor2; % nJ/nm
spectrum2_588 = spectrum2_588*1e-9/(Nt*dt*1e-12); % W/nm
spectrum_588 = spectrum1_588 + spectrum2_588;
spectrum_588 = spectrum_588./max(spectrum_588);
figure;
plot(lambda,spectrum_588,'linewidth',6,'Color','k');
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
set(gca,'fontsize',40);
xlim([1000,1060]);
print(gcf,'f_Qswitching.pdf','-dpdf');
% noise
spectrum1_681 = abs(fftshift(ifft(output_field{1}(:,:,681)),1)).^2*factor_correct_unit1; % nJ/THz
spectrum1_681 = spectrum1_681.*factor1*correction_factor1; % nJ/nm
spectrum1_681 = spectrum1_681*1e-9/(gain_rate_eqn.woPulse.Nt*dt_woPulse*1e-12); % W/nm
spectrum1_681 = interp1(lambda1,spectrum1_681,lambda,'linear',0);
spectrum2_681 = abs(fftshift(ifft(output_field{2}(:,:,681)),1)).^2*factor_correct_unit2; % nJ/THz
spectrum2_681 = spectrum2_681.*factor2; % nJ/nm
spectrum2_681 = spectrum2_681*1e-9/(Nt*dt*1e-12); % W/nm
spectrum_681 = spectrum1_681 + spectrum2_681;
spectrum_681 = spectrum_681./max(spectrum_681);
figure;
plot(lambda,spectrum_681,'linewidth',6,'Color','k');
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
set(gca,'fontsize',40);
xlim([1000,1060]);
print(gcf,'f_noise.pdf','-dpdf');

load(sprintf('oscillator_%u.mat',4000));
spectrum1_001 = abs(fftshift(ifft(output_field{1}(:,:,001)),1)).^2*factor_correct_unit1; % nJ/THz
spectrum1_001 = spectrum1_001.*factor1*correction_factor1; % nJ/nm
spectrum1_001 = spectrum1_001*1e-9/(gain_rate_eqn.woPulse.Nt*dt_woPulse*1e-12); % W/nm
spectrum1_001 = interp1(lambda1,spectrum1_001,lambda,'linear',0);
spectrum2_001 = abs(fftshift(ifft(output_field{2}(:,:,001)),1)).^2*factor_correct_unit2; % nJ/THz
spectrum2_001 = spectrum2_001.*factor2; % nJ/nm
spectrum2_001 = spectrum2_001*1e-9/(Nt*dt*1e-12); % W/nm
spectrum_001 = spectrum1_001 + spectrum2_001;
spectrum_001 = spectrum_001./max(spectrum_001);
figure;
plot(lambda,spectrum_001,'linewidth',6,'Color','k');
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
set(gca,'fontsize',40);
xlim([1000,1060]);
print(gcf,'f_modelocking.pdf','-dpdf');

%% Q-switching burst
load(sprintf('oscillator_%u.mat',1000));

idx = 588 + (-20:20);
t_Q = t1 + (0:length(idx)-1)*(prop_output_woPulse.dt*gain_rate_eqn.woPulse.Nt); t_Q = t_Q(:); t_Q = t_Q - t_Q(1);
Q_power = feval(@(x) x(:), abs(output_field{1}(:,:,idx)).^2);

figure;
plot(t_Q/1e6,Q_power,'linewidth',2,'Color','b');
xlim([0,max(t_Q)/1e6]);
xlabel('Time (\mus)');
ylabel('Power (W)');
set(gca,'fontsize',20);
print(gcf,'Q-switched_burst.pdf','-dpdf');