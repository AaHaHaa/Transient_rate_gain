clearvars; close all;

tr_data = load('transient_ASE_evolutions.mat');
ss_data = load('./With steady-state rate-gain code/steady-state_ASE_evolutions.mat');

%% Power at each z
tr_ASE = struct('forward', squeeze(mean(abs(tr_data.prop_output(1).fields.forward).^2,1))*1e3,... % mW
                'backward',squeeze(mean(abs(tr_data.prop_output(1).fields.backward).^2,1))*1e3);  % mW

ss_ASE = struct('forward', squeeze(trapz(fftshift(ss_data.f,1),ss_data.prop_output.Power.ASE.forward))*1e3,...
                'backward',squeeze(trapz(fftshift(ss_data.f,1),ss_data.prop_output.Power.ASE.backward))*1e3);

figure;
h = plot(tr_data.prop_output(1).z,[tr_ASE.forward,tr_ASE.backward],'linewidth',3);
set(h(1),'Color','b'); set(h(2),'Color','r');
hold on;
h2 = plot(ss_data.prop_output(1).z,[ss_ASE.forward,ss_ASE.backward],'MarkerSize',7,'Marker','o','MarkerEdgeColor','k','LineStyle','None');
set(h2(1),'MarkerFaceColor','b'); set(h2(2),'MarkerFaceColor','r');
hold off;
set(gca,'fontsize',20);
xlim([tr_data.prop_output(1).z(1),tr_data.prop_output(1).z(end)]);
xlabel('Propagation (m)'); ylabel('Power (mW)');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);
l = legend('forward (tr)','backward (tr)','forward (ss)','backward (ss)');
set(gca,'Position',[0.085,0.187,0.88,0.73])
print(gcf,'ASE power comparison.pdf','-dpdf','-bestfit');

%% Temporal profiles
tr_Nt1 = length(tr_data.t);
tr_dt1 = tr_data.prop_output(1).dt;
tr_t1 = (-tr_Nt1/2:tr_Nt1/2-1)'*tr_dt1; % ps
tr_t1 = tr_t1 - min(tr_t1);

figure;
tr_ASE_forward = abs(tr_data.prop_output(1).fields.forward(:,:,end)).^2*1e3;
avg_tr_ASE_forward = smooth(tr_ASE_forward,30);
h = plot(tr_t1/1e6,[tr_ASE_forward,avg_tr_ASE_forward,tr_ASE.forward(end)*ones(tr_Nt1,1)],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','w','linewidth',4,'linestyle','--');
set(gca,'fontsize',20);
xlabel('Time (ms)'); ylabel('Power (mW)');
legend('|A(t)|^2','smoothed |A(t)|^2','ss value');
print(gcf,'tr_ASE_forward(t).pdf','-dpdf');
print(gcf,'tr_ASE_forward(t).jpg','-djpeg');

figure;
tr_ASE_backward = abs(tr_data.prop_output(1).fields.backward(:,:,1)).^2*1e3;
avg_tr_ASE_backward = smooth(tr_ASE_backward,30);
h = plot(tr_t1/1e6,[tr_ASE_backward,avg_tr_ASE_backward,tr_ASE.backward(1)*ones(tr_Nt1,1)],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','w','linewidth',4,'linestyle','--');
set(gca,'fontsize',20);
xlabel('Time (ms)'); ylabel('Power (mW)');
legend('|A(t)|^2','smoothed |A(t)|^2','ss value');
print(gcf,'tr_ASE_backward(t).pdf','-dpdf');
print(gcf,'tr_ASE_backward(t).jpg','-djpeg');

%% Spectra
% steady-state
c = 299792.458; % nm/ps;
ss_lambda = c./fftshift(ss_data.f,1); % nm
ss_factor = c./ss_lambda.^2; % change the spectrum from frequency domain into wavelength domain
ss_ASEf_forward = ss_data.prop_output.Power.ASE.forward(:,:,end).*ss_factor; % W/nm
ss_ASEf_backward = ss_data.prop_output.Power.ASE.backward(:,:,1).*ss_factor; % W/nm

% transient
tr_f1_real = tr_data.f; % THz
tr_lambda1 = c./tr_f1_real; % nm
tr_factor_correct_unit = (tr_Nt1*tr_dt1)^2/1e3; % nJ/THz
tr_factor = c./tr_lambda1.^2; % change the spectrum from frequency domain into wavelength domain
tr_f1 = (-tr_Nt1/2:tr_Nt1/2-1)'/(tr_Nt1*tr_dt1); % THz
cs = tr_data.prop_output(1).dt/tr_data.sim.incoherent_dt;

tr_ASEf_forward = abs(fftshift(ifft(tr_data.prop_output(1).fields.forward(:,:,end)),1)).^2*tr_factor_correct_unit; % nJ/THz
tr_ASEf_forward = tr_ASEf_forward.*tr_factor; % nJ/nm
tr_ASEf_forward = tr_ASEf_forward*1e-9/(tr_Nt1*tr_dt1*1e-12)/cs; % W/nm
tr_ASEf_backward = abs(fftshift(ifft(tr_data.prop_output(1).fields.backward(:,:,1)),1)).^2*tr_factor_correct_unit; % nJ/THz
tr_ASEf_backward = tr_ASEf_backward.*tr_factor; % nJ/nm
tr_ASEf_backward = tr_ASEf_backward*1e-9/(tr_Nt1*tr_dt1*1e-12)/cs; % W/nm

figure;
plot(tr_lambda1,tr_ASEf_forward*1e3,'linewidth',2,'Color','b');
hold on;
plot(ss_lambda,ss_ASEf_forward*1e3,'MarkerSize',7,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','None');
hold off;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)'); ylabel('PSD (mW/nm)');
legend('forward (tr)','forward (ss)');
xlim([1000,1120]);
print('ASE_forward(lambda).pdf','-dpdf');
print('ASE_forward(lambda).jpg','-djpeg');

figure;
plot(tr_lambda1,tr_ASEf_backward*1e3,'linewidth',2,'Color','b');
hold on;
plot(ss_lambda,ss_ASEf_backward*1e3,'MarkerSize',7,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','None');
hold off;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)'); ylabel('PSD (mW/nm)');
legend('backward (tr)','backward (ss)');
xlim([1000,1120]);
print('ASE_backward(lambda).pdf','-dpdf');
print('ASE_backward(lambda).jpg','-djpeg');

% spectra at window 2 (with the pulse)
tr_Nt2 = length(tr_data.f);
tr_dt2 = mean(diff(tr_data.t));
tr_f2 = tr_data.f;
tr_lambda2 = c./tr_f2; % nm
tr_factor_correct_unit = (tr_Nt2*tr_dt2)^2/1e3; % nJ/THz
tr_factor = c./tr_lambda2.^2; % change the spectrum from frequency domain into wavelength domain
tr_ASEf_forward2 = abs(fftshift(ifft(tr_data.prop_output(2).fields.forward(:,:,end)),1)).^2*tr_factor_correct_unit; % nJ/THz
tr_ASEf_forward2 = tr_ASEf_forward2.*tr_factor; % nJ/nm
tr_ASEf_forward2 = tr_ASEf_forward2*1e-9/(tr_Nt2*tr_dt2*1e-12); % W/nm
figure;
plot(tr_lambda2,tr_ASEf_forward2*1e3,'linewidth',2,'Color','b');
hold on;
plot(ss_lambda,ss_ASEf_forward*1e3,'MarkerSize',7,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','None');
hold off;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)'); ylabel('PSD (mW/nm)');
legend('forward (tr)','forward (ss)');
xlim([1000,1120]);

tr_ASEf_backward2 = abs(fftshift(ifft(tr_data.prop_output(2).fields.backward(:,:,1)),1)).^2*tr_factor_correct_unit; % nJ/THz
tr_ASEf_backward2 = tr_ASEf_backward2.*tr_factor; % nJ/nm
tr_ASEf_backward2 = tr_ASEf_backward2*1e-9/(tr_Nt2*tr_dt2*1e-12); % W/nm
figure;
plot(tr_lambda2,tr_ASEf_backward2*1e3,'linewidth',2,'Color','b');
hold on;
plot(ss_lambda,ss_ASEf_backward*1e3,'MarkerSize',7,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','LineStyle','None');
hold off;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)'); ylabel('PSD (mW/nm)');
legend('backward (tr)','backward (ss)');
xlim([1000,1120]);