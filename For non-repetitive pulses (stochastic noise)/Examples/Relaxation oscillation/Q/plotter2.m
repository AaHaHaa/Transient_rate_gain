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

%% Q-switching burst
load(sprintf('oscillator_%u.mat',6000));

idx =267 + (-50:50);
t_Q = t1 + (0:length(idx)-1)*(prop_output_woPulse.dt*gain_rate_eqn.woPulse.Nt); t_Q = t_Q(:); t_Q = t_Q - t_Q(1);
Q_power = feval(@(x) x(:), abs(output_field{1}(:,:,idx)).^2);

figure;
plot(t_Q/1e6,Q_power,'linewidth',1,'Color','b');
xlim([0,max(t_Q)/1e6]);
xlabel('Time (\mus)');
ylabel('Power (W)');
set(gca,'fontsize',20);
print(gcf,'Q-switched_burst.pdf','-dpdf');

% closeview
figure;
plot(t_Q/1e6,Q_power,'linewidth',2,'Color','b');
xlim([0.39,0.45]);
set(gca,'fontsize',30);
print(gcf,'Q-switched_burst (closeview).pdf','-dpdf');