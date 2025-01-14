function [At_out,Power_out,...
          T_delay_out,...
          N] = SteppingCaller_rategain(sim,gain_rate_eqn,...
                                       num_zPoints,save_points,num_zPoints_persave,...
                                       initial_condition,...
                                       extended_n2_prefactor,...
                                       SK_info,SRa_info, SRb_info,...
                                       extended_expD,...
                                       extended_haw, extended_hbw,...
                                       At_noise)
%STEPPINGCALLER_RATEGAIN It attains the field after propagation inside the  
%gain medium solved by the rate equations.
%   
% The computation of this code is based on
%   1. Lindberg et al., "Accurate modeling of high-repetition rate ultrashort pulse amplification in optical fibers", Scientific Reports (2016)
%   2. Chen et al., "Optimization of femtosecond Yb-doped fiber amplifiers for high-quality pulse compression", Opt. Experss (2012)
%   3. Gong et al., "Numerical modeling of transverse mode competition in strongly pumped multimode fiber lasers and amplifiers", Opt. Express (2007)
%
%   Please go to "gain_info.m" file to find the information about some input arguments.
%   The information of some other input arguments are inside "Transient_gain_UPPE_propagate.m"
%

num_windows = length(initial_condition);
Nt = size(initial_condition(1).fields.forward,1);
dt = [initial_condition.dt];

t = zeros(Nt*num_windows,1);
t_each_window = cell(1,num_windows);
cumsum_Nt = Nt*(1:num_windows);
t_each_window{1} = (0:Nt-1)'*dt(1);
t(1:Nt) = t_each_window{1};
for window_j = 2:num_windows
    t_each_window{window_j} = t(cumsum_Nt(window_j-1))+dt(window_j-1) + (0:Nt-1)'*dt(window_j);
    t(cumsum_Nt(window_j-1) + (1:Nt)) = t_each_window{window_j};
end

N_guess = cell(1,1,1,num_windows);

%% Pump direction
if gain_rate_eqn.copump_power == 0
    if gain_rate_eqn.counterpump_power == 0
        gain_rate_eqn.pump_direction = 'co'; % use 'co' for zero pump power
    else
        gain_rate_eqn.pump_direction = 'counter';
    end
else
    if gain_rate_eqn.counterpump_power == 0
        gain_rate_eqn.pump_direction = 'co';
    else
        gain_rate_eqn.pump_direction = 'bi';
    end
end

%% Segment the entire propagation
% In situations where iterations are needed, such as counterpumping, all
% the information is saved during pulse propagation.
% However, if GPU is used, this might overload its memory. In such a
% situation, the propagation is segmented into several pieces. Only one
% piece is computed by GPU at each time instance while the rest is kept
% into RAM.

if ~sim.gpu_yes % use CPU
    num_segments = 1;
    segments = num_zPoints; % the number of z points of all segments; [num_segment1,num_segment2...]
    zPoints_each_segment = num_zPoints;
else % use GPU
    % Below calculates how to segment the data based on the predicted size (memory usage) of the total data.
    precision = 8; % "double" precision
    mem_complex_number = precision*2;
    % The size of the variables:
    variable_size.A = num_windows*(Nt*num_windows)*10*num_zPoints; % forward and backward
    variable_size.Power_pump = num_windows*(Nt*num_windows)*num_zPoints; % forward and backward
    variable_size.cross_sections = 10*(Nt*num_windows);
    variable_size.overlap_factor = numel(gain_rate_eqn.overlap_factor);
    variable_size.N_total = numel(gain_rate_eqn.N_total);
    variable_size.N = (Nt*num_windows)*(length(gain_rate_eqn.energy_levels)-1)*num_zPoints;
    var_field = fieldnames(variable_size);
    used_memory = 0;
    for i = 1:length(var_field)
        used_memory = used_memory + variable_size.(var_field{i});
    end
    used_memory = used_memory*mem_complex_number;

    num_segments = ceil(used_memory/gain_rate_eqn.memory_limit);
    if num_segments == 1
        segments = num_zPoints; % the number of z points of all segments; [num_segment1,num_segment2...]
    else
        zPoints_each_segment = ceil(num_zPoints/num_segments);
        num_segments = ceil(num_zPoints/zPoints_each_segment); % this operation is super important because a small increase of zPoints_each_segment (due to the ceiling function) can reduce num_segments by a few
        segments = [zPoints_each_segment*ones(1,num_segments-1) num_zPoints-zPoints_each_segment*(num_segments-1)]; % the number of z points of all segments; [num_segment1,num_segment2...]
        if segments(end) == 0
            segments = segments(1:end-1);
        end
    end
end

%% Initialization
% They're both in the frequency domain.
if isequal(gain_rate_eqn.pump_direction,'co')
    segment_idx = 1;
else % 'counter', 'bi'      
    segment_idx = num_segments;
end

% Initialize everything in the specified segment, including both propagating forward and backward
% Everything is initialized because it's not sure which will be used.
% If forward propagation is done first, we need the backward data since in this code, both directions are, by default, considered.
% If backward propagation is done first, we need the forward data since in this code, both directions are, by default, considered.
[At_forward,At_backward,...
 Power_pump_forward,Power_pump_backward] = initialization('both',...
                                                          gain_rate_eqn,...
                                                          segment_idx,num_segments,...
                                                          Nt,num_windows,...
                                                          num_zPoints,...
                                                          initial_condition);

%% Propagations
% =========================================================================
% Start the first pulse propagation
% =========================================================================
% If counter/bi-pumping, backward-propagate first without the signal to set up the population-inversion level.
if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
                  gain_propagate('backward',At_noise);
end

% Start the forward pulse propagation
[T_delay_out,N] = gain_propagate('forward',At_noise);

% -------------------------------------------------------------------------
% Initialize some variables for the following iterations
% -------------------------------------------------------------------------
if sim.gpu_yes
    pulse_energy = zeros(1,gain_rate_eqn.max_iterations,'gpuArray');
    
    if gain_rate_eqn.verbose
        ASE_forward = zeros(1,gain_rate_eqn.max_iterations,'gpuArray');
        ASE_backward = zeros(1,gain_rate_eqn.max_iterations,'gpuArray');
    end
else
    pulse_energy = zeros(1,gain_rate_eqn.max_iterations);
    
    if gain_rate_eqn.verbose
        ASE_forward = zeros(1,gain_rate_eqn.max_iterations);
        ASE_backward = zeros(1,gain_rate_eqn.max_iterations);
    end
end

% =========================================================================
% Start the iterations (the 2nd, 3rd, ... pulse propagations) if necessary
% =========================================================================
% There's no need to compute backward propagation if copumping and ignoring ASE,
% so the previous first forward propagation solves everything.
if ~isequal(gain_rate_eqn.pump_direction,'co') || gain_rate_eqn.include_ASE
    pulse_energy(1) = calc_total_energy(At_forward{1,1,end,2},dt(2)); % nJ
    if gain_rate_eqn.verbose
        fprintf('Gain rate equation, iteration %u: pulse energy = %7.6g(nJ)\n',1,pulse_energy(1));
        
        % ASE powers when there is no pulse
        ASE_forward(1) = mean(abs(At_forward{1,1,end,1}).^2,1)*1e3; % mW
        ASE_backward(1) = mean(abs(At_backward{1,1,1,1}).^2,1)*1e3; % mW
        fprintf('                    iteration %u: forward  ASE power (at seed output end) = %7.6g(mW)\n',1,ASE_forward(1));
        fprintf('                    iteration %u: backward ASE power (at seed input end)  = %7.6g(mW)\n',1,ASE_backward(1));
    end
    finish_iteration = false;
    for i = 2:gain_rate_eqn.max_iterations
        % -----------------------------------------------------------------
        % Propagate back and forth
        % -----------------------------------------------------------------
        % Backward propagation
                          gain_propagate('backward',At_noise);

        % Forward propagation
        [T_delay_out,N] = gain_propagate('forward',At_noise);

        % -----------------------------------------------------------------
        % Check convergence
        % -----------------------------------------------------------------
        pulse_energy(i) = calc_total_energy(At_forward{1,1,end,2},dt(2)); % nJ
        if gain_rate_eqn.verbose
            fprintf('Gain rate equation, iteration %u: pulse energy = %7.6g(nJ)\n',i,pulse_energy(i));
            
            % ASE powers when there is no pulse
            ASE_forward(i) = mean(abs(At_forward{1,1,end,1}).^2,1)*1e3; % mW
            ASE_backward(i) = mean(abs(At_backward{1,1,1,1}).^2,1)*1e3; % mW
            fprintf('                    iteration %u: forward  ASE power (at seed output end) = %7.6g(mW)\n',i,ASE_forward(i));
            fprintf('                    iteration %u: backward ASE power (at seed input end)  = %7.6g(mW)\n',i,ASE_backward(i));
        end
        if pulse_energy(i-1)==0 || ... % no propagating pulse; this is to check the validity of the following pulse_energy convergence check due to the 1/pulse_energy(i-1) factor
           abs((pulse_energy(i)-pulse_energy(i-1))./pulse_energy(i-1)) < gain_rate_eqn.tol % pulse reaches a steady state
            finish_iteration = true;
        end
        % Plot the convergence
        if gain_rate_eqn.verbose && i > 10 && ~(i==10+1 && finish_iteration) % this final condition means that the iteration finishes right when it's the 11th iteration. Don't plot anything in this situation.
            if i == 10+1 % plot it the first time; initialize the figure
                fig_gain_iterations = figure('Name','Gain iterations');
            end
            figure(fig_gain_iterations);
            plot(1:i,pulse_energy(1:i),'linewidth',2,'Color','b');
            set(gca,'YColor','b');
            xlabel('Iterations'); ylabel('Pulse energy (nJ)');
            set(gca,'fontsize',20);
            drawnow;
            % Close the convergence plot when it's done
            if finish_iteration
                close(fig_gain_iterations);
            end
        end
        
        % -----------------------------------------------------------------
        % Done!
        % -----------------------------------------------------------------
        if finish_iteration
            break;
        end
        % -----------------------------------------------------------------
        % Error! Iterations of the gain-fiber computation isn't finished
        % within "max_iterations" iterations.
        % -----------------------------------------------------------------
        if i == gain_rate_eqn.max_iterations
            warning('UPPE_rategain:NotConvergedError',...
                    ['The iteration of forward and backward propagation of the gain fiber doesn''t converge to a steady state within %u iterations.\n',...
                     'Please run it again with a larger number of maximum iterations.'],gain_rate_eqn.max_iterations);
            break;
        end
    end
end

%% Output:
saved_zPoints = 1:num_zPoints_persave:num_zPoints;

% Transform them into arrays
At_out = cell(1,num_windows);
Power_out = cell(1,num_windows);
for window_j = 1:num_windows
    At_out{window_j}   = struct('forward', cell2mat(At_forward (1,1,saved_zPoints,window_j)),...
                                'backward',cell2mat(At_backward(1,1,saved_zPoints,window_j)));

    Power_out{window_j} = struct('pump',struct('forward', cell2mat(Power_pump_forward (1,1,saved_zPoints,window_j)),...
                                               'backward',cell2mat(Power_pump_backward(1,1,saved_zPoints,window_j))));
end

N = shiftdim(N,2);

%%
function [T_delay_out,N] = gain_propagate(direction,At_noise)
%GAIN_PROPAGATE Runs the corresponding propagation method based on "direction".

if isempty(N_guess{1})
    switch direction
        case 'forward'
            At_forward_guess = At_forward(:,:,1,:);
            At_backward_guess = At_backward(:,:,1,:);
            pump_power = gain_rate_eqn.copump_power;
        case 'backward'
            At_forward_guess = At_forward(:,:,end,:);
            At_backward_guess = At_backward(:,:,end,:);
            pump_power = gain_rate_eqn.counterpump_power;
    end
    
    N_guess_value = solve_gain_rate_eqn_steadyState(gain_rate_eqn,sim.gpu_yes,...
                                                    At_forward_guess,At_backward_guess,...
                                                    pump_power,...
                                                    Nt,dt);
    for window_i = 1:num_windows
        if sim.gpu_yes
            N_guess{window_i} = ones(Nt,1,'gpuArray').*N_guess_value.'; % initial guess for solving the population during propagation
        else
            N_guess{window_i} = ones(Nt,1).*N_guess_value.'; % initial guess for solving the population during propagation
        end
    end
end

T_delay_out = zeros(save_points,1);

zPoints_each_segment = segments(1);

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
    
    % Because circshift is slow on GPU, I discard it.
    %last_result = ifft(circshift(initial_condition.fields,-tCenter),[],1);
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

% Initialization
switch direction
    case 'forward'
        [At_forward,Power_pump_forward] = initialization('forward',...
                                                         gain_rate_eqn,...
                                                         1,num_segments,...
                                                         Nt,num_windows,...
                                                         num_zPoints,...
                                                         initial_condition);
        % Put the first segment into GPU
        if sim.gpu_yes
            [At_forward,At_backward,...
             Power_pump_forward,Power_pump_backward] = mygpuArray2(1,segments(1),...
                                                                   At_forward,At_backward,...
                                                                   Power_pump_forward,Power_pump_backward);
        end
    case 'backward'
        [At_backward,Power_pump_backward] = initialization('backward',...
                                                           gain_rate_eqn,...
                                                           num_segments,num_segments,...
                                                           Nt,num_windows,...
                                                           num_zPoints,...
                                                           initial_condition);
        % Put the last segment into GPU
        if sim.gpu_yes
            [At_forward,At_backward,...
             Power_pump_forward,Power_pump_backward] = mygpuArray2(num_zPoints-segments(end)+1,num_zPoints,...
                                                                   At_forward,At_backward,...
                                                                   Power_pump_forward,Power_pump_backward);
        end
end

% Initialize N to be exported, the ion density of the upper state
N = cell(1,1,1,num_windows);
for window_i = 1:num_windows
    N{window_i} = zeros(Nt,length(gain_rate_eqn.energy_levels)-1,save_points);
end

if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('Transient_gain_UPPE_propagate:ProgressBarNameError',...
              '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.0%%',sim.progress_bar_name),...
        'Name',sprintf('Running UPPE (%s): %s...',direction,sim.progress_bar_name),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h_progress_bar,'canceling',0);
    
    % Create the cleanup object
    cleanupObj = onCleanup(@()cleanMeUp(h_progress_bar));
    
    % Use this to control the number of updated time for the progress bar below 1000 times.
    count_progress_bar = 1;
    num_progress_updates = 1000;
end

% Then start the propagation
last_N = N_guess;
dN = zeros(length(gain_rate_eqn.energy_levels)-1,1);
for ii = 2:num_zPoints
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('Transient_gain_UPPE_propagate:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end
    
    % =====================================================================
    % UPPE: Run the correct step function depending on the options chosen.
    % =====================================================================
    switch direction
        % -----------------------------------------------------------------
        % Forward propagation
        % -----------------------------------------------------------------
        case 'forward'
            Zi = ii; % z-index of the forward-propagating distance
            
            % Load the initial powers and field
            if Zi == 2 % the first/starting z-index
                last_At_forward = At_forward(:,:,1,:); % = initial_condition(1:2).fields.forward
                last_Power_pump_forward = Power_pump_forward(:,:,1,:);
            end
            last_At_backward = At_backward(:,:,Zi-1,:); % = initial_condition(1:2).fields.forward
            last_Power_pump_backward = Power_pump_backward(:,:,Zi-1,:);
            
            % Since backward pump power and the backward field are a function of time,
            % temporal offset from pulse-centering should be considered for them as well.
            if sim.pulse_centering && (~isnan(TCenter) && TCenter ~= 0)
                last_Power_pump_backward_allWindows = cat(1,last_Power_pump_backward{:});
                
                for window_i = 1:num_windows
                    last_Power_pump_backward{window_i} = interp1(t,last_Power_pump_backward_allWindows,t_each_window{window_i}+TCenter*dt(2),'linear','extrap');
                    if sim.gpu_yes && mod(window_i,2) == 1 % pseudo-field window doesn't have a lot of points, so CPU is faster
                        last_Power_pump_backward{window_i} = gather(last_Power_pump_backward{window_i});
                    end

                    if mod(window_i,2) == 0
                        if TCenter > 0
                            last_At_backward{window_i} = [last_At_backward{window_i}(1+TCenter:end,:);last_At_backward{window_i}(1:TCenter,:)];
                        elseif TCenter < 0
                            last_At_backward{window_i} = [last_At_backward{window_i}(end+1+TCenter:end,:);last_At_backward{window_i}(1:end+TCenter,:)];
                        end
                    end
                end
            end
            
            [last_At_forward,...
             last_Power_pump_forward,...
             last_N,dN] = forward_stepping_RK4IP_rategain(last_At_forward,last_At_backward,...
                                                          last_Power_pump_forward,last_Power_pump_backward,...
                                                          last_N,dN,...
                                                          Nt,dt,num_windows,...
                                                          sim,gain_rate_eqn,...
                                                          SK_info,SRa_info,SRb_info,...
                                                          extended_haw,extended_hbw,...
                                                          extended_n2_prefactor,extended_expD,...
                                                          At_noise(1,:));
            
            % Update "forward" only
            At_forward(:,:,Zi,:) = last_At_forward;
            Power_pump_forward(:,:,Zi,:) = last_Power_pump_forward;
            
            % Save N
            if Zi == num_zPoints % save the last N
                last_N = find_periodic_transient_N('forward',...
                                                   sim,gain_rate_eqn,...
                                                   Nt,dt,num_windows,...
                                                   last_N,dN,...
                                                   last_At_forward,last_At_backward,...
                                                   last_Power_pump_forward,last_Power_pump_backward);

                if sim.gpu_yes
                    for window_i = 1:num_windows
                        N{window_i}(:,:,end) = gather(last_N{window_i}); % (Nt,num_levels-1)
                    end
                else
                    for window_i = 1:num_windows
                        N{window_i}(:,:,end) = last_N{window_i}; % (Nt,num_levels-1)
                    end
                end

                for window_i = 1:num_windows
                    N{window_i} = N{window_i}/gather(gain_rate_eqn.N_total);
                end
            end
            % Current N is computed by the next propagating step, so there is Zi-2, not Zi-1.
            % Coincidentally, this command also helps save the first N.
            if rem(Zi-2, num_zPoints_persave) == 0
                if sim.gpu_yes
                    for window_i = 1:num_windows
                        N{window_i}(:,:,int64((Zi-2)/num_zPoints_persave+1)) = gather(last_N{window_i});
                    end
                else
                    for window_i = 1:num_windows
                        N{window_i}(:,:,int64((Zi-2)/num_zPoints_persave+1)) = last_N{window_i};
                    end
                end
            end
            
            % Whenever at the end of each segment, gather the current data
            % from GPU back to the RAM and then throw those in the next
            % segment to GPU.
            if rem(Zi,zPoints_each_segment) == 0 || Zi == num_zPoints
                ready_for_next_segment_into_GPU = true;
            else
                ready_for_next_segment_into_GPU = false;
            end
        % -----------------------------------------------------------------
        % Backward propagation
        % -----------------------------------------------------------------
        case 'backward'
            Zi = num_zPoints+1 - ii; % z-index of the backward-propagating distance
            
            % Load the final powers
            if Zi == num_zPoints+1 - 2 % ii=2; the first/starting z-index
                last_Power_pump_backward = Power_pump_backward(:,:,num_zPoints,:);
                last_At_backward = At_backward(:,:,num_zPoints,:);
            end
            last_Power_pump_forward = Power_pump_forward(:,:,Zi+1,:);
            last_At_forward = At_forward(:,:,Zi+1,:);
            
            [last_At_backward,...
             last_Power_pump_backward,...
             last_N,dN] = backward_stepping_RK4IP_rategain(last_At_forward,last_At_backward,...
                                                           last_Power_pump_forward,last_Power_pump_backward,...
                                                           last_N,dN,...
                                                           Nt,dt,num_windows,...
                                                           sim,gain_rate_eqn,...
                                                           At_noise(2,:));
            
            % Update "backward" only
            At_backward(:,:,Zi,:) = last_At_backward;
            Power_pump_backward(:,:,Zi,:) = last_Power_pump_backward;
            
            % Whenever at the beginning of each segment, gather the current
            % data from GPU back to the RAM and then throw those in the
            % next segment to GPU.
            if rem(Zi,zPoints_each_segment) == 1
                ready_for_next_segment_into_GPU = true;
            else
                ready_for_next_segment_into_GPU = false;
            end
    end
    
    % =====================================================================
    % Some post-stepping checks, saves, and updates
    % =====================================================================
    if isequal(direction,'forward')
        % Check for any NaN elements
        for window_i = 1:num_windows
            if any(any(isnan(last_At_forward{window_i})))
                error('SteppingCaller_rategain:NaNError',...
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
                temporal_profile = abs(last_At_forward{window_i}).^2;
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
                        if TCenter > 0
                            last_At_forward{window_i} = [last_At_forward{window_i}(1+TCenter:end,:);last_At_forward{window_i}(1:TCenter,:)];
                            At_noise{1,window_i} = cat(1,At_noise{1,window_i}(1+TCenter:end,:,:),At_noise{1,window_i}(1:TCenter,:,:));
                        elseif TCenter < 0
                            last_At_forward{window_i} = [last_At_forward{window_i}(end+1+TCenter:end,:);last_At_forward{window_i}(1:end+TCenter,:)];
                            At_noise{1,window_i} = cat(1,At_noise{1,window_i}(end+1+TCenter:end,:,:),At_noise{1,window_i}(1:end+TCenter,:,:));
                        end
                    end
                end

                if sim.gpu_yes
                    T_delay = T_delay + gather(TCenter*dt(2));
                else
                    T_delay = T_delay + TCenter*dt(2);
                end
            end
            if rem(Zi-1, num_zPoints_persave) == 0
                T_delay_out(int64((Zi-1)/num_zPoints_persave)+1) = T_delay;
            end
        end
    end
    
    % Gather the current GPU data back to RAM and those in the next segment from RAM to GPU
    if ready_for_next_segment_into_GPU
        current_segment_idx = ceil(Zi/zPoints_each_segment);
        switch direction
            case 'forward'
                next_segment_idx = current_segment_idx + 1;
            case 'backward'
                next_segment_idx = current_segment_idx - 1;
        end

        cumsum_segments = [0,cumsum(segments)];
        start_idx = cumsum_segments(current_segment_idx)+1;
        end_idx   = cumsum_segments(current_segment_idx+1);
        [At_forward,At_backward,...
         Power_pump_forward,Power_pump_backward] = mygather2(start_idx,end_idx,...
                                                             At_forward,At_backward,...
                                                             Power_pump_forward,Power_pump_backward);

        if next_segment_idx > 0 && next_segment_idx <= length(segments)
            start_idx = cumsum_segments(next_segment_idx)+1;
            end_idx   = cumsum_segments(next_segment_idx+1);
            [At_forward,At_backward,...
             Power_pump_forward,Power_pump_backward] = mygpuArray2(start_idx,end_idx,...
                                                                   At_forward,At_backward,...
                                                                   Power_pump_forward,Power_pump_backward);
        end
    end
    
    % Report current status in the progress bar's message field
    if sim.progress_bar
        if num_zPoints < num_progress_updates || floor((ii-1)/((num_zPoints-1)/num_progress_updates)) == count_progress_bar
            waitbar((ii-1)/(num_zPoints-1),h_progress_bar,sprintf('%s%6.1f%%',sim.progress_bar_name,(ii-1)/(num_zPoints-1)*100));
            count_progress_bar = count_progress_bar+1;
        end
    end
end

N_guess = last_N;

end

end

%% initialization
function varargout = initialization(direction,...
                                    gain_rate_eqn,...
                                    segment_idx,num_segments,...
                                    Nt,num_windows,...
                                    zPoints,...
                                    initial_condition)
%INITIALIZATION initializes "signal_fields" and "Powers" based on
%"segment_idx/num_segment".
%
%   If segment_idx = 1, they need to include copump_power, initial_fields, and initial forward ASE.
%   If segment_idx = num_segment(the last segment), "Power" needs to include counterpump_power if it's nonzero.
%                                                   They also need to include initial backward ASE.
%   Otherwise, they're just structures with zero matrices.
%
%   The reason I use cell arrays instead of a matrix (N,num_modes,zPoints) for signal_fields and Power:
%       It's faster!
%       e.g. "signal_fields(:,:,zi) = signal_fields_next" is very slow.

    function output = initialize_zeros()
        output = cell(1,1,zPoints,num_windows);
        output(:) = {zeros(Nt,1)};
    end

% =========================================================================
% Initialization
% =========================================================================
if ismember(direction,{'forward','both'})
    Power_pump_forward = initialize_zeros();
    At_forward = initialize_zeros();
end
if ismember(direction,{'backward','both'})
    Power_pump_backward = initialize_zeros();
    At_backward = initialize_zeros();
end

% -------------------------------------------------------------------------
% Put in the necessary information
% -------------------------------------------------------------------------
% Pump power
if ismember(direction,{'forward','both'})
    if segment_idx == 1
        for window_i = 1:num_windows
            At_forward{1,1,1,window_i} = initial_condition(window_i).fields.forward;
            
            if ismember(gain_rate_eqn.pump_direction,{'co','bi'})
                Power_pump_forward{1,1,1,window_i} = gain_rate_eqn.copump_power*ones(Nt,1);
            end
        end
    end
end
if ismember(direction,{'backward','both'})
    if segment_idx == num_segments
        for window_i = 1:num_windows
            At_backward{1,1,end,window_i} = initial_condition(window_i).fields.backward;

            if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
                Power_pump_backward{1,1,end,window_i} = gain_rate_eqn.counterpump_power*ones(Nt,1);
            end
        end
    end
end

% =========================================================================
% Output arguments
% =========================================================================
switch direction
    case 'forward'
        varargout = {At_forward,Power_pump_forward};
    case 'backward'
        varargout = {At_backward,Power_pump_backward};
    otherwise % 'both'
        varargout = {At_forward,At_backward,...
                     Power_pump_forward,Power_pump_backward};
end

end

%% CALC_TOTAL_ENERGY
function total_energy = calc_total_energy(At,dt)
%CALC_TOTAL_ENERGY

total_energy = sum(abs(At).^2)*dt/1e3; % nJ;

end

%% MYGPUARRAY2
function varargout = mygpuArray2(starti,endi,varargin)
%MYGPUARRAY2 It throws a part of the inputs to GPU, from "starti" to "endi"
% Input arguments should be in "cell" arrays.

varargout = varargin;
for i = 1:nargin-2
    varargout{i}(:,:,starti:endi,2:2:end) = cellfun(@gpuArray,varargin{i}(:,:,starti:endi,2:2:end),'UniformOutput',false);
end

end

%% MYGATHER2
function varargout = mygather2(starti,endi,varargin)
%MYGATHER2 It gathers a part of the inputs from GPU to the RAM, from "starti" to "endi"
% Input arguments should be in "cell" arrays.

varargout = varargin;
for i = 1:nargin-2
    varargout{i}(:,:,starti:endi,2:2:end) = cellfun(@gather,varargin{i}(:,:,starti:endi,2:2:end),'UniformOutput',false);
end

end

%% CLEANMEUP
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end