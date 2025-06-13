close all; clearvars;

%rmpath('../../../../../GMMNLSE/GMMNLSE algorithm/','../../../../../GMMNLSE/user_helpers/');
%rmpath('../../../../For repetitive pulses (stochastic noise)/UPPE algorithm/','../../../../For repetitive pulses (stochastic noise)/user_helpers/');
addpath('../../../UPPE algorithm/','../../../user_helpers/');

%%
load('../steady state/ss_for_initial_N.mat');

%%
sim.cuda_dir_path = '../cuda';
sim.dz = sim.save_period;
sim.progress_bar = false;
%sim.gpu_yes = true;
sim.cs = 1;

fiber.n2 = 2.3e-20;

%%
sim.incoherent.dt = dt/3;
incoherent_lambda0 = 1070e-9; % m
sim.incoherent.f0 = 2.99792458e-4/incoherent_lambda0;
incoherent_f = sim.incoherent.f0+(-Nt/2:Nt/2-1)'/(Nt*sim.incoherent.dt); % THz
incoherent_lambda = c./(incoherent_f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( sim,gain_rate_eqn,{ifftshift(incoherent_lambda,1),ifftshift(lambda,1)} );

%% Setup initial conditions
previous_population = permute(prop_output.population(end,:,:,:),[1,4,3,2]);

rep_rate = 1/gain_rate_eqn.t_rep; % Hz

tfwhm = 1; % ps
total_energy = 1e-10; % nJ

prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

prop_output.fields = struct('forward',prop_output.fields);
prop_output.population = previous_population;

tfwhm_woPulse = 10e3;
dt_woPulse = (1/rep_rate*1e12-time_window)/Nt;
time_window_woPulse = dt_woPulse*Nt;
total_energy_woPulse = 1e-10;

pulse_lambda0 = 1030e-9;
f1 = incoherent_f; % THz
f1 = ifftshift(linspace(f1(2),f1(1),Nt)',1); % THz
f_now = f1(1);
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = (f_pulse - sim.incoherent.f0)/(dt_woPulse/sim.incoherent.dt);
prop_output_woPulse = build_noisy_MMgaussian(tfwhm_woPulse, tfwhm_woPulse*1, time_window_woPulse, total_energy_woPulse, total_energy_woPulse, 1, Nt, 0.005, {'ifft',freq_shift});
prop_output_woPulse.fields = struct('forward',prop_output_woPulse.fields);
prop_output_woPulse.population = previous_population;

initial_condition(1) = prop_output_woPulse;
initial_condition(2) = prop_output;

%%
t_woPulse = (-Nt/2:Nt/2-1)'*prop_output_woPulse.dt; % ps

%% cavity parameters
max_rt = 5000; % maximum roundtrips (in case it doesn't converge)
each_max_rt = 1000;
OC = 0.8; % output coupling
loss = 0.5; % the total loss of the cavity
%saturation_power = 2500; % the saturation power of the saturable absorber; W
moddepth = 1; % the modulation depth of the saturable absorber
tol_convergence = 1e-3; % the tolerance of the convergence of the output ANDi pulse energy

% Spectral filter parameters
gaussexpo = 1;
plot_filter_yes = false;%true;
spectral_filter = struct('bw',6, ... % bandwidth (nm)
                         'cw',1030); % center wavelength (nm)

%% Save the general information
save('oscillator_general_info.mat');

%% Run the simulation
initial_condition(2).fields.forward = fft(ifft(initial_condition(2).fields.forward(:,:,end)).*exp(1i*0.2/2*2*pi*ifftshift(f-sim.f0).^2));
prop_output = Transient_gain_UPPE_propagate(fiber, initial_condition, sim, gain_rate_eqn);

% -----------------------------------------------------------------
% Output coupler
prop_output(1).fields.forward = sqrt(1-OC)*prop_output(1).fields.forward(:,:,end);
prop_output(2).fields.forward = sqrt(1-OC)*prop_output(2).fields.forward(:,:,end);

% -----------------------------------------------------------------
% Saturable absorber
tmp1 = struct('dt',prop_output(1).dt,'fields',prop_output(1).fields.forward);
tmp2 = struct('dt',prop_output(2).dt,'fields',prop_output(2).fields.forward);
saturation_power_Q = 0.0005;%min(2500,max(abs(tmp2.fields).^2));
tmp1 = saturable_absorber_action_simple(tmp1, saturation_power_Q, moddepth);
saturation_power = min(2500,max(abs(tmp2.fields).^2));
tmp2 = saturable_absorber_action_simple(tmp2, saturation_power, moddepth);

% -----------------------------------------------------------------
% Spectral filter
tmp11 = tmp1; tmp11.dt = sim.incoherent.dt;
tmp1 = gaussian_spectral_filter(tmp1, sim.incoherent.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter
tmp2 = gaussian_spectral_filter(tmp2, sim.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter

tmp1.fields = sqrt(1-loss)*tmp1.fields(:,:,end);
tmp2.fields = sqrt(1-loss)*tmp2.fields(:,:,end);

prop_output(1).fields.forward = tmp1.fields;
prop_output(2).fields.forward = tmp2.fields;
rt_num = 1;
while rt_num <= max_rt
    output_field = {zeros(Nt,1,each_max_rt),zeros(Nt,1,each_max_rt)}; % the output pulse
    output_energy = zeros(each_max_rt,2);
    population1 = zeros(each_max_rt,1);
    population2 = zeros(each_max_rt,1);
    population3 = zeros(each_max_rt,1);
    population4 = zeros(each_max_rt,1);
    pump1 = zeros(each_max_rt,1);
    pump2 = zeros(each_max_rt,1);
    pump3 = zeros(each_max_rt,1);
    pump4 = zeros(each_max_rt,1);


    rt_num_thousand = mod(rt_num,each_max_rt);

    while rt_num_thousand ~= 0
        rt_num_thousand = mod(rt_num,each_max_rt);
        
        t_iteration_start = tic;
        fprintf('Iteration %d\n', rt_num);
        % -----------------------------------------------------------------
        previous_population = prop_output(2).population(end,:,:);
        last_fields1 = prop_output(1).fields.forward(:,:,end);
        last_fields2 = prop_output(2).fields.forward(:,:,end);
        dt1 = prop_output(1).dt;
        dt2 = prop_output(2).dt;
        clearvars prop_output;
        prop_output(1).population = previous_population;
        prop_output(2).population = zeros(size(previous_population));
        prop_output(1).fields.forward = last_fields1;
        prop_output(2).fields.forward = last_fields2;
        prop_output(1).dt = dt1;
        prop_output(2).dt = dt2;
    
        prop_output(2).fields.forward = fft(ifft(prop_output(2).fields.forward(:,:,end)).*exp(1i*0.2/2*2*pi*ifftshift(f-sim.f0).^2));
        prop_output = Transient_gain_UPPE_propagate(fiber, prop_output, sim, gain_rate_eqn);
    
        if rt_num_thousand == 0
            idx = each_max_rt;
        else
            idx = rt_num_thousand;
        end
        population1(idx) = prop_output(2).population(end,:,floor(fiber.L0/sim.save_period/4)+1);
        population2(idx) = prop_output(2).population(end,:,floor(fiber.L0/sim.save_period/4*2)+1);
        population3(idx) = prop_output(2).population(end,:,floor(fiber.L0/sim.save_period/4*3)+1);
        population4(idx) = prop_output(2).population(end,:,floor(fiber.L0/sim.save_period/4*4)+1);

        pump1(idx) = prop_output(2).Power.pump.forward(end,:,floor(fiber.L0/sim.save_period/4)+1);
        pump2(idx) = prop_output(2).Power.pump.forward(end,:,floor(fiber.L0/sim.save_period/4*2)+1);
        pump3(idx) = prop_output(2).Power.pump.forward(end,:,floor(fiber.L0/sim.save_period/4*3)+1);
        pump4(idx) = prop_output(2).Power.pump.forward(end,:,floor(fiber.L0/sim.save_period/4*4)+1);
    
        if mod(rt_num,100)==0
            close all;
        end

        %{
        figure;
        yyaxis left;
        plot(t1/1e6,prop_output(1).population(:,:,end)*100,'linewidth',2);
        ylabel('N_1 (%)');
        yyaxis right;
        plot(t1/1e6,abs(prop_output(1).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
        xlabel('Time (\mus)'); ylabel('Power (kW)');
        set(gca,'fontsize',20);
        
        figure;
        yyaxis left;
        plot(t,prop_output(2).population(:,:,end)*100,'linewidth',2);
        ylabel('N_1 (%)');
        yyaxis right;
        plot(t,abs(prop_output(2).fields.forward(:,:,end)).^2/1e3,'linewidth',2);
        xlabel('Time (ps)'); ylabel('Power (kW)');
        set(gca,'fontsize',20);
    
        drawnow;
        %}
        % -----------------------------------------------------------------
        % Output coupler
        output_field{1}(:,:,idx) = sqrt(OC)*prop_output(1).fields.forward(:,:,end);
        output_field{2}(:,:,idx) = sqrt(OC)*prop_output(2).fields.forward(:,:,end);
        prop_output(1).fields.forward = sqrt(1-OC)*prop_output(1).fields.forward(:,:,end);
        prop_output(2).fields.forward = sqrt(1-OC)*prop_output(2).fields.forward(:,:,end);
        
        % Energy of the output field
        output_energy(idx,1) = trapz(abs(output_field{1}(:,:,idx)).^2).*prop_output(1).dt/1e3; % energy in nJ
        output_energy(idx,2) = trapz(abs(output_field{2}(:,:,idx)).^2).*prop_output(2).dt/1e3;
    
        if mod(rt_num,100)==0
            figure;
            yyaxis left;
            plot(output_energy(1:idx,1));
            yyaxis right;
            plot(output_energy(1:idx,2));
            title('Energy')
        
            figure;
            yyaxis left;
            plot(population1(1:idx)*100);
            yyaxis right;
            plot(population3(1:idx)*100);
            title('N_1');
    
            figure;
            plot(pump1(1:idx)); hold on;
            plot(pump2(1:idx));
            plot(pump3(1:idx));
            plot(pump4(1:idx)); hold off;
            title('Pump')

            drawnow;
        end
    
        % -----------------------------------------------------------------
        % Saturable absorber
        tmp1 = struct('dt',prop_output(1).dt,'fields',prop_output(1).fields.forward);
        tmp2 = struct('dt',prop_output(2).dt,'fields',prop_output(2).fields.forward);
        if rt_num < 50
            saturation_power_Q = 0.0001;
        else
            saturation_power_Q = 1;
        end
        tmp1 = saturable_absorber_action_simple(tmp1, saturation_power_Q, moddepth);
        saturation_power = min(2500,max(abs(tmp2.fields).^2));
        tmp2 = saturable_absorber_action_simple(tmp2, saturation_power, moddepth);
        
        % -----------------------------------------------------------------
        % Spectral filter
        tmp11 = tmp1; tmp11.dt = sim.incoherent.dt;
        tmp1 = gaussian_spectral_filter(tmp1, sim.incoherent.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter
        tmp2 = gaussian_spectral_filter(tmp2, sim.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter
        
        tmp1.fields = sqrt(1-loss)*tmp1.fields(:,:,end);
        tmp2.fields = sqrt(1-loss)*tmp2.fields(:,:,end);
    
        prop_output(1).fields.forward = tmp1.fields;
        prop_output(2).fields.forward = tmp2.fields;
        
        % -----------------------------------------------------------------
    
        rt_num = rt_num + 1;
    end

    close all;
    save(sprintf('oscillator_%u.mat',rt_num-1),'output_field','output_energy','population1','population2','population3','population4','pump1','pump2','pump3','pump4');
end