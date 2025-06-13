clearvars; close all;

%% Paper data
Fig2_data = readmatrix('one_nsPulse (Fig 2).csv');
Fig3_data = readmatrix('one_nsPulse (Fig 3).csv');
Fig5_data = readmatrix('two_nsPulse_offsetFreq (Fig 5).csv');
Fig6_data = readmatrix('two_nsPulse_offsetTime (Fig 6).csv');
Fig7_data = readmatrix('one_fsPulse (Fig 7).csv');

%% Simulation data
Fig23_sim = load('one_nsPulse.mat');
Fig5_sim = load('two_nsPulses_offsetFreq.mat');
Fig6_sim = load('two_nsPulses_offsetTime.mat');
Fig7_sim = load('one_fsPulse.mat');

%% Plot
% Figs. 2 and 3
Fig2_P = abs(Fig23_sim.prop_output(2).fields.forward(:,:,end)).^2;
figure;
plot(Fig23_sim.t/1e3,Fig2_P/1e3,'linewidth',2,'Color','b');
hold on;
plot(Fig2_data(:,1),Fig2_data(:,2),'linewidth',2,'Color','r');
hold off;
set(gca,'fontsize',20);
xlim([-1,1]);
xlabel('Time (ns)');
ylabel('Power (kW)');
legend('Ours','Prior work');
print(gcf,'Fig2.pdf','-dpdf');

factor = 1./Fig23_sim.lambda.^2; % change the spectrum from frequency domain into wavelength domain
Fig3_spectrum = abs(fftshift(ifft(Fig23_sim.prop_output(2).fields.forward(:,:,end)),1)).^2.*factor; Fig3_spectrum = Fig3_spectrum/max(Fig3_spectrum);
figure;
plot(Fig23_sim.lambda,Fig3_spectrum,'linewidth',2,'Color','b');
hold on;
plot(Fig3_data(:,1),Fig3_data(:,2)/max(Fig3_data(:,2)),'linewidth',2,'Color','r');
hold off;
set(gca,'fontsize',20);
xlim([1020,1040]);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
l = legend('Ours','Prior work'); set(l,'location','northwest');
print(gcf,'Fig3.pdf','-dpdf');

% Fig. 5
factor = 1./Fig5_sim.lambda.^2; % change the spectrum from frequency domain into wavelength domain
Fig5_spectrum = abs(fftshift(ifft(Fig5_sim.prop_output(2).fields.forward(:,:,end)),1)).^2.*factor; Fig5_spectrum = Fig5_spectrum/max(Fig5_spectrum);
figure;
plot(Fig5_sim.lambda,Fig5_spectrum,'linewidth',2,'Color','b');
hold on;
plot(Fig5_data(:,1),Fig5_data(:,2)/max(Fig5_data(:,2)),'linewidth',2,'Color','r');
hold off;
set(gca,'fontsize',20);
xlim([1010,1060]);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
legend('Ours','Prior work');
print(gcf,'Fig5.pdf','-dpdf');

% Fig. 6
Fig6_P = circshift(abs(Fig6_sim.prop_output(2).fields.forward(:,:,end)).^2,ceil(Fig6_sim.prop_output(2).t_delay(end)/Fig6_sim.prop_output(2).dt));
figure;
plot(Fig6_sim.t/1e3,Fig6_P/1e3,'linewidth',2,'Color','b');
hold on;
plot(Fig6_data(:,1),Fig6_data(:,2),'linewidth',2,'Color','r');
hold off;
xlim([-1,1]);
set(gca,'fontsize',20);
xlabel('Time (ns)');
ylabel('Power (kW)');
legend('Ours','Prior work');
print(gcf,'Fig6.pdf','-dpdf');

% Fig. 7
factor = 1./Fig7_sim.lambda.^2; % change the spectrum from frequency domain into wavelength domain
Fig7_spectrum = abs(fftshift(ifft(Fig7_sim.prop_output(2).fields.forward(:,:,end)),1)).^2.*factor; Fig7_spectrum = Fig7_spectrum/max(Fig7_spectrum);
figure;
plot(Fig7_sim.lambda,Fig7_spectrum,'linewidth',2,'Color','b');
hold on;
plot(Fig7_data(:,1),Fig7_data(:,2)/max(Fig7_data(:,2)),'linewidth',2,'Color','r');
hold off;
set(gca,'fontsize',20);
xlim([1000,1090]);
xlabel('Wavelength (nm)');
ylabel('PSD (norm.)');
legend('Ours','Prior work');
print(gcf,'Fig7.pdf','-dpdf');