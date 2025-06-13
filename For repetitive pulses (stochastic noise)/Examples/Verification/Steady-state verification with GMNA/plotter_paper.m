clearvars; close all;

tr_data = load('transient_GMNA.mat');
ss_data = load('./With steady-state rate-gain code/steady-state_GMNA.mat');

tr_data.Nt = length(tr_data.t);
tr_data.dt = mean(diff(tr_data.t));
ss_data.Nt = length(ss_data.t);
ss_data.dt = mean(diff(ss_data.t));

%% Pulse energy and Power at each z
tr_energy = squeeze(trapz(tr_data.t,abs(tr_data.prop_output(2).fields.forward).^2,1))/1e3; % nJ
ss_energy = squeeze(trapz(ss_data.t,abs(ss_data.prop_output.fields).^2,1))/1e3; % nJ

figure;
h = plot(tr_data.prop_output(1).z,[tr_energy,ss_energy]);
set(h(1),'Color','b','linewidth',3); set(h(2),'MarkerFaceColor','r','MarkerSize',7,'Marker','o','MarkerEdgeColor','k','LineStyle','None');
set(gca,'fontsize',20);
xlim([tr_data.prop_output(1).z(1),tr_data.prop_output(1).z(end)]);
xlabel('Propagation (m)'); ylabel('Pulse energy (nJ)');
l = legend('tr','ss'); set(l,'location','northwest');
print(gcf,'pulse energy comparison.pdf','-dpdf');

%% Population
tr_N1 = squeeze(tr_data.prop_output(2).population(1,:,:,:));
ss_N1 = squeeze(ss_data.prop_output.population(1,:,:,:));

figure;
h = plot(tr_data.prop_output(1).z,[tr_N1,ss_N1]*100);
set(h(1),'Color','b','linewidth',3); set(h(2),'MarkerFaceColor','r','MarkerSize',7,'Marker','o','MarkerEdgeColor','k','LineStyle','None');
set(gca,'fontsize',20);
xlim([tr_data.prop_output(1).z(1),tr_data.prop_output(1).z(end)]);
xlabel('Propagation (m)'); ylabel('N_1 (%)');
legend('tr','ss');
print(gcf,'N1 comparison.pdf','-dpdf');

%% Temporal profiles
downsampling = 30;
figure;
h(1) = plot(tr_data.t,abs(tr_data.prop_output(2).fields.forward(:,:,end)).^2/1e3);
hold on;
h(2) = plot(ss_data.t(1:downsampling:end),abs(ss_data.prop_output.fields(1:downsampling:end,:,end)).^2/1e3);
hold off;
set(h(1),'Color','b','linewidth',3); set(h(2),'MarkerFaceColor','r','MarkerSize',7,'Marker','o','MarkerEdgeColor','k','LineStyle','None');
set(gca,'fontsize',20);
xlim([-5,5]);
xlabel('Time (ps)'); ylabel('Power (kW)');
legend('tr','ss');
print(gcf,'A2 comparison.pdf','-dpdf');

%% Spectra
% steady-state
c = 299792.458; % nm/ps;
ss_lambda = c./ss_data.f; % nm
ss_factor_correct_unit = (ss_data.Nt*ss_data.dt)^2/1e3; % nJ/THz
ss_factor = c./ss_lambda.^2; % change the spectrum from frequency domain into wavelength domain
ss_Af_forward = abs(fftshift(ifft(ss_data.prop_output.fields(:,:,end)),1)).^2*ss_factor_correct_unit.*ss_factor; % nJ/nm

% transient
tr_lambda = c./tr_data.f; % nm
tr_factor_correct_unit = (tr_data.Nt*tr_data.dt)^2/1e3; % nJ/THz
tr_factor = c./tr_lambda.^2; % change the spectrum from frequency domain into wavelength domain
tr_Af_forward = abs(fftshift(ifft(tr_data.prop_output(2).fields.forward(:,:,end)),1)).^2*tr_factor_correct_unit.*tr_factor; % nJ/nm

figure;
plot(tr_lambda,tr_Af_forward,'linewidth',2,'Color','b');
hold on;
plot(ss_lambda(1:downsampling*3:end),ss_Af_forward(1:downsampling*3:end),'MarkerSize',7,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','None');
hold off;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)'); ylabel('PSD (nJ/nm)');
legend('tr','ss');
xlim([980,1200]);
print('A(lambda).pdf','-dpdf');