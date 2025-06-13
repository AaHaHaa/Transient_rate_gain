function [At_out,Power_out,...
          save_z,save_dz,...
          T_delay_out,...
          N] = SteppingCaller_adaptive_rategain(sim,gain_rate_eqn,...
                                                save_z0,save_points,...
                                                initial_condition,...
                                                extended_n2_prefactor,...
                                                extended_D_op,...
                                                SK_info,SRa_info,SRb_info,...
                                                extended_haw,extended_hbw,...
                                                At_noise)
%STEPPINGCALLER_ADAPTIVE_RATEGAIN It attains the field after propagation 
%inside the gain medium solved by the rate equations with the adaptive-step
%RK4IP algorithm.
%   
%   Please go to "gain_info.m" file to find the information about some input arguments.
%   The info of some other input arguments are inside "Transient_gain_UPPE_propagate.m"
%

gain_rate_eqn.pump_direction = 'co';

%% Propagations
[At_forward,Power_pump_forward,...
 save_z,save_dz,...
 T_delay_out,...
 N] = gain_propagate(sim,gain_rate_eqn,...
                     save_points,save_z0,...
                     extended_n2_prefactor,extended_D_op,...
                     SK_info,SRa_info,SRb_info,extended_haw,extended_hbw,...
                     At_noise,...
                     initial_condition);

%% Output:
if sim.gpu_yes
    [At_forward,Power_pump_forward] = mygather(At_forward,Power_pump_forward);
end

% This is for consistency with the non-adaptive-step computation which has 
% At_out{1}.forward, At_out{1}.backward, At_out{2}.forward, and At_out{2}.backward.
% Here, in the adaptive-step computation, there is only At_out{2}.forward.
num_windows = length(initial_condition);
At_out = cell(1,num_windows);
Power_out = cell(1,num_windows);
for window_i = 2:2:num_windows
    At_out{window_i-1} = struct('forward',[]); % adaptive-step method doesn't consider ASE, so there is nothing in this window
    At_out{window_i}   = struct('forward',cell2mat(At_forward(1,1,:,window_i))); % window with coherent field

    Power_out{window_i-1}.pump.forward = cell2mat(Power_pump_forward(1,1,:,window_i-1));
    Power_out{window_i  }.pump.forward = cell2mat(Power_pump_forward(1,1,:,window_i  ));
end

N = shiftdim(N,2);

end

%%
function [signal_fields,Power_pump_forward,...
          save_z,save_dz,...
          T_delay_out,...
          N] = gain_propagate(sim,gain_rate_eqn,...
                              save_points,save_z,...
                              extended_n2_prefactor,extended_D_op,...
                              SK_info,SRa_info,SRb_info,...
                              extended_haw,extended_hbw,...
                              At_noise,...
                              initial_condition)
%GAIN_PROPAGATE Runs the corresponding propagation method based on "direction".

num_windows = length(initial_condition);
Nt = size(initial_condition(2).fields.forward,1);
dt = [initial_condition.dt];

t = zeros(sum(Nt),1);
t_each_window = cell(1,num_windows);
cumsum_Nt = Nt*(1:num_windows);
t_each_window{1} = (0:Nt-1)'*dt(1);
t(1:Nt) = t_each_window{1};
for window_i = 2:num_windows
    t_each_window{window_i} = t(cumsum_Nt(window_i-1))+dt(window_i-1) + (0:Nt-1)'*dt(window_i);
    t(cumsum_Nt(window_i-1) + (1:Nt)) = t_each_window{window_i};
end

dummy_var = cell(1,1,1,num_windows);
for window_i = 1:num_windows
    if sim.gpu_yes
        dummy_var{window_i} = zeros(Nt,1,'gpuArray');
    else
        dummy_var{window_i} = zeros(Nt,1);
    end
end

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);

% Pulse centering based on the moment of its intensity
if sim.pulse_centering
    % Center the pulses
    TCenter_each_window = zeros(1,num_windows/2);
    for window_i = 2:2:num_windows
        temporal_profile = abs(initial_condition(window_i).fields.forward).^2;
        temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
        TCenter_each_window(window_i/2) = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,1)/sum(temporal_profile,1));
    end
    TCenter = floor(mean(TCenter_each_window));

    if ~isnan(TCenter) && TCenter ~= 0
        for window_i = 2:2:num_windows
            if TCenter > 0
                initial_condition(window_i).fields.forward = [initial_condition(window_i).fields.forward(1+TCenter:end,:);initial_condition(window_i).fields.forward(1:TCenter,:)];
            elseif TCenter < 0
                initial_condition(window_i).fields.forward = [initial_condition(window_i).fields.forward(end+1+TCenter:end,:);initial_condition(window_i).fields.forward(1:end+TCenter,:)];
            end
        end

        if sim.gpu_yes
            TCenter = gather(TCenter);
        end
        T_delay = TCenter*dt(2);
    else
        T_delay = 0;
    end
else
    T_delay = 0;
end
T_delay_out(1) = T_delay;

[signal_fields,Power_pump_forward] = initialization(sim,gain_rate_eqn,...
                                                    Nt,num_windows,save_points,...
                                                    initial_condition);
last_Power_pump_forward = Power_pump_forward(:,:,1,:);
last_signal_fields = signal_fields(:,:,1,:); % = initial_condition(2:2:end).fields.forward

% Initialize N to be exported, the ion density of the upper state
N = cell(1,1,1,num_windows);
for window_i = 1:num_windows
    N{window_i} = zeros(Nt,length(gain_rate_eqn.energy_levels)-1,save_points);
end

if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('SteppingCaller_adaptive_rategain:ProgressBarNameError',...
              '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.0%%',sim.progress_bar_name),...
        'Name',sprintf('Running UPPE with transient gain: %s...',sim.progress_bar_name),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h_progress_bar,'canceling',0);
    
    % Create the cleanup object
    cleanupObj = onCleanup(@()cleanMeUp(h_progress_bar));
    
    % Use this to control the number of updated time for the progress bar below 1000 times.
    num_progress_updates = 1000;
    progress_bar_z = (1:num_progress_updates)*save_z(end)/num_progress_updates;
    progress_bar_i = 1;
end

% max dz
if ~isfield(sim.adaptive_dz,'max_dz')
    sim.adaptive_dz.max_dz = sim.save_period/10;
end

sim.dz = min(1e-6,sim.adaptive_dz.max_dz); % m; start with a small value to avoid initial blowup
save_dz(1) = sim.dz;

% Then start the propagation
z = 0;
save_i = 2; % the 1st one is the initial field
a5 = cell(1,1,1,num_windows); % the temporary variable in the forward propagation

last_N = cell(1,1,1,num_windows);
N_guess = solve_gain_rate_eqn_steadyState(gain_rate_eqn,sim.gpu_yes,...
                                          last_signal_fields,[],...
                                          gain_rate_eqn.copump_power,...
                                          Nt,dt);
for window_i = 1:num_windows
    if sim.gpu_yes && mod(window_i,2) == 0
        last_N{window_i} = ones(Nt,1,'gpuArray').*N_guess.'; % initial guess for solving the population during propagation
    else
        last_N{window_i} = ones(Nt,1).*N_guess.'; % initial guess for solving the population during propagation
    end
end
dN = zeros(length(gain_rate_eqn.energy_levels)-1,1);
sim.last_dz = 1; % randomly put a number, 1, for initialization
while z+eps(z) < save_z(end) % eps(z) here is necessary due to the numerical error
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('SteppingCaller_adaptive_rategain:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end
    
    ever_fail = false;
    previous_signal_fields = last_signal_fields;
    previous_a5 = a5;

    success = false;
    while ~success
        if ever_fail
            last_signal_fields = previous_signal_fields;
            a5 = previous_a5;
        end

        [last_signal_fields,a5,...
         last_Power_pump_forward,...
         last_N,dN,...
         opt_dz,success] = stepping_RK4IP_rategain_adaptive(last_signal_fields,last_Power_pump_forward,...
                                                       last_N,dN,...
                                                       Nt,dt,num_windows,...
                                                       sim,gain_rate_eqn,...
                                                       SK_info,SRa_info,SRb_info,...
                                                       extended_haw,extended_hbw,...
                                                       extended_n2_prefactor,extended_D_op,...
                                                       At_noise,...
                                                       a5,...
                                                       dummy_var);

        if ~success
            if opt_dz < 1e-10
                error('SteppingCaller_adaptive_rategain:adaptiveRK4IPError',...
                      'Adaptive RK4IP continues to fail.\nCheck simulation parameters.');
            end

            ever_fail = true;

            sim.dz = opt_dz;
        end
    end
    sim.last_dz = sim.dz; % previous dz

    % Check for any NaN elements
    for window_i = 2:2:num_windows
        if any(any(isnan(last_signal_fields{window_i/2})))
            error('SteppingCaller_adaptive_rategain:NaNError',...
                  'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
        end
    end

    % Center the pulse in the time window
    % Important:
    % In the modified shot-noise approach, the noise cannot be changed, so it needs to be translated as well.
    % This took me more than 12 hours of debugging to realize it.
    % Otherwise, the output field, if it has a strong frequency shift and shifts a lot in time relative to the time window, 
    % the noise without translation cannot match with the translated field,
    % resulting in a different noise field overlapped with the coherent pulse.
    % This will artificially create a noisy output.
    if sim.pulse_centering
        TCenter_each_window = zeros(1,num_windows/2);
        for window_i = 2:2:num_windows
            temporal_profile = abs(last_signal_fields{window_i}).^2;
            temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
            TCenter_each_window(window_i/2) = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,1)/sum(temporal_profile,1));
        end
        TCenter = floor(mean(TCenter_each_window));
        if ~isnan(TCenter) && TCenter ~= 0
            last_Power_pump_forward_allWindows = cat(1,last_Power_pump_forward{:});
            last_N_allWindows = cat(1,last_N{:});
            
            for window_i = 1:num_windows
                last_Power_pump_forward{window_i} = interp1(t,last_Power_pump_forward_allWindows,t_each_window{window_i}+TCenter*dt(2),'linear','extrap');
                last_N{window_i} = interp1(t,last_N_allWindows,t_each_window{window_i}+TCenter*dt(2),'linear','extrap');
                if sim.gpu_yes && mod(window_i,2) == 1
                    last_Power_pump_forward{window_i} = gather(last_Power_pump_forward{window_i});
                    last_N{window_i} = gather(last_N{window_i});
                end

                if mod(window_i,2) == 0
                    a5{window_i} = fft(a5{window_i},[],1);
                
                    if TCenter > 0
                        a5{window_i} = ifft([a5{window_i}(1+TCenter:end,:);a5{window_i}(1:TCenter,:)],[],1);
                        last_signal_fields{window_i} = [last_signal_fields{window_i}(1+TCenter:end,:);last_signal_fields{window_i}(1:TCenter,:)];
                        At_noise{window_i} = cat(1,At_noise{window_i}(1+TCenter:end,:,:),At_noise{window_i}(1:TCenter,:,:));
                    elseif TCenter < 0
                        a5{window_i} = ifft([a5{window_i}(end+1+TCenter:end,:);a5{window_i}(1:end+TCenter,:)],[],1);
                        last_signal_fields{window_i} = [last_signal_fields{window_i}(end+1+TCenter:end,:);last_signal_fields{window_i}(1:end+TCenter,:)];
                        At_noise{window_i} = cat(1,At_noise{window_i}(end+1+TCenter:end,:,:),At_noise{window_i}(1:end+TCenter,:,:));
                    end
                end
            end
            
            if sim.gpu_yes
                TCenter = gather(TCenter);
            end
            T_delay = T_delay + TCenter*dt(2);
        end
    end

    % remember the old dz for scaling dN to fast find the optimal dN such
    % that N+dN satisfies the periodic boundary condition for the next
    % stepping computation at z+dz. Since dz is small, dN should scale
    % linearly with dz. Because this is an adaptive-step computation,
    % scaling of dN is necessary.
    old_dz = sim.dz;

    % Update z
    z = z + sim.dz;
    % Because the adaptive-step algorithm determines the step size by 
    % checking the error of the spectral intensities from RK3 and RK4, it
    % might ignore changes at the weakest part of the spectrum. This
    % happens in cases of noise-seeded four-wave mixing and noise-seeded
    % Raman scattering far from the pulse central frequency.
    %
    % To account for this effect, I limit the step size to be 10x smaller 
    % than the effective maximum beat length which is
    % 2*pi/(max(eff_betas)-min(eff_betas)).
    % eff_betas is from betas, propagation constants, throughout the time 
    % window but scaled according to the spectral intensity to prevent 
    % taking into account betas where there's no light at all or where 
    % there's some light starting to grow.
    eff_range_D = zeros(1,num_windows/2);
    for window_i = 2:2:num_windows
        extended_last_signal_fields_w = ifft(last_signal_fields{window_i},gain_rate_eqn.acyclic_conv_stretch(Nt),1);
        eff_range_D(window_i/2) = find_range_D(abs(extended_last_signal_fields_w).^2,imag(extended_D_op));
    end
    min_beat_length = 2*pi/max(eff_range_D);
    dz_resolve_beat_length = min_beat_length/4;
    sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz,dz_resolve_beat_length]);

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if z == sim.last_dz
        ready_save_N = true;
    end
    if z >= save_z(end)-eps(z) % save the last N
        last_N = find_periodic_transient_N('forward',...
                                           sim,gain_rate_eqn,...
                                           Nt,dt,num_windows,...
                                           last_N,dN,...
                                           last_signal_fields,dummy_var,...
                                           last_Power_pump_forward,dummy_var);

        if sim.gpu_yes
            for window_i = 1:num_windows
                N{window_i}(:,:,end) = gather(last_N{window_i}); % (Nt,num_levels-1)
            end
        else
            for window_i = 1:num_windows
                N{window_i}(:,:,end) = last_N{window_i}; % (Nt,num_levels-1)
            end
        end

        N = cellfun(@(x)x/gather(gain_rate_eqn.N_total),N,'UniformOutput',false);
    end
    % Below saves N, the population of each energy level:
    % Current N is computed by the next propagating step, so 
    % there is "ready_save_N" is to delay the saving process 
    % by one propagating step with the "save_i-1" command.
    % Coincidentally, this command also helps save the first
    % N.
    if ready_save_N
        if sim.gpu_yes
            for window_i = 1:num_windows
                N{window_i}(:,:,save_i-1) = gather(last_N{window_i}); % (Nt,num_levels-1)
            end
        else
            for window_i = 1:num_windows
                N{window_i}(:,:,save_i-1) = last_N{window_i}; % (Nt,num_levels-1)
            end
        end
        ready_save_N = false; % finish saving, then reset it back to false
    end
    if z >= save_z(save_i)-eps(z)
        if sim.gpu_yes
            save_z(save_i) = gather(z);
            save_dz(save_i) = gather(sim.last_dz);
            Power_pump_forward(1,1,save_i,:) = cellfun(@gather,last_Power_pump_forward,'UniformOutput',false);
            signal_fields(1,1,save_i,:) = cellfun(@gather,last_signal_fields,'UniformOutput',false);
        else
            save_z(save_i) = z;
            save_dz(save_i) = sim.last_dz;
            Power_pump_forward(1,1,save_i,:) = last_Power_pump_forward;
            signal_fields(1,1,save_i,:) = last_signal_fields;
        end
        ready_save_N = true;

        T_delay_out(save_i) = T_delay;

        save_i = save_i + 1;
    end
    
    % Scale dN according to varying dz
    dN = dN*sim.dz/old_dz;

    % Report current status in the progress bar's message field
    if sim.progress_bar
        if z >= progress_bar_z(progress_bar_i)
            waitbar(gather(z/save_z(end)),h_progress_bar,sprintf('%s%6.1f%%',sim.progress_bar_name,z/save_z(end)*100));
            progress_bar_i = find(z<progress_bar_z,1);
        end
    end
end

end

%% Initialization
function [At_forward,Power_pump_forward] = initialization(sim,gain_rate_eqn,...
                                                          Nt,num_windows,save_points,...
                                                          initial_condition)
%INITIALIZATION initializes "signal_fields" and "Power"
%
%   The reason I use cell arrays instead of a matrix (N,num_modes,z_points) for signal_fields and Power:
%       It's faster!
%       e.g. "signal_fields(:,:,zi) = signal_fields_next" is very slow.

    function output = initialize_zeros()
        output = cell(1,1,save_points,num_windows);
        output(:) = {zeros(Nt,1)};
    end

% =========================================================================
% Initialization
% =========================================================================
Power_pump_forward  = initialize_zeros();
At_forward = initialize_zeros();

% -------------------------------------------------------------------------
% Put in the necessary information
% -------------------------------------------------------------------------
for window_i = 1:num_windows
    Power_pump_forward{1,1,1,window_i} = gain_rate_eqn.copump_power*ones(Nt,1);
end
% Signal field
for window_i = 2:2:num_windows
    At_forward{1,1,1,window_i} = initial_condition(window_i).fields.forward; % with coherent pulses
end

% -------------------------------------------------------------------------
% GPU
if sim.gpu_yes
    [At_forward,Power_pump_forward] = mygpuArray(At_forward,Power_pump_forward);
end

end

%% MYGPUARRAY
function varargout = mygpuArray(varargin)
%MYGPUARRAY

varargout =varargin;
for i = 1:nargin
    varargout{i}(:,:,:,2:2:end) = cellfun(@gpuArray,varargin{i}(:,:,:,2:2:end),'UniformOutput',false);
end

end

%% MYGATHER
function varargout = mygather(varargin)
%MYGATHER

varargout = cell(1,nargin);
for i = 1:nargin
    varargout{i} = cellfun(@gather,varargin{i},'UniformOutput',false);
end

end

%% EFF_RANGE_D
function eff_range_D = find_range_D(spectrum,D)
%FIND_RANGE_D
%
% For an adaptive-dz method, the maximum dz is also limited by the 
% range of the propagation constant, beta0.
% If the FWM, Raman, or anything else happens for multiple frequencies 
% where dz can't resolve their beta0 difference, the outcome can be 
% wrong.
% Here, I multiply the (intensity)^(1/2) of the spectrum to the beta0 to 
% consider beta0 difference of the pulse and exclude those without the 
% pulse but within the frequency window. (1/2) is to maximize the 
% contribution of the weak part of the spectrum.

nonzero_field = max(spectrum)~=0;
spectrum = spectrum./max(spectrum(:));

eff_D = D(:,nonzero_field).*(spectrum(:,nonzero_field)).^(1/2); % I use ^(1/2) to emphasize the weak part
eff_range_D = max(eff_D(:)) - min(eff_D(:));

end

%% CLEANMEUP
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end