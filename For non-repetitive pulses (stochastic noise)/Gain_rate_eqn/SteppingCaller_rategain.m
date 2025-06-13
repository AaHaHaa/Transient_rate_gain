function [At_out,Power_out,...
          T_delay,...
          N_out] = SteppingCaller_rategain(sim,gain_rate_eqn,...
                                           num_zPoints,...
                                           initial_condition,...
                                           extended_n2_prefactor,...
                                           SK_info,SRa_info, SRb_info,...
                                           extended_expD,...
                                           extended_haw, extended_hbw,...
                                           At_noise)
%STEPPINGCALLER_RATEGAIN It attains the field after propagation inside the  
%gain medium solved by the rate equations.
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

%% Pump direction
gain_rate_eqn.pump_direction = 'co';

%% Segment the entire propagation
% In situations where iterations are needed, such as counterpumping, all
% the information is saved during pulse propagation.
% However, if GPU is used, this might overload its memory. In such a
% situation, the propagation is segmented into several pieces. Only one
% piece is computed by GPU at each time instance while the rest is kept
% into RAM.

if ~sim.gpu_yes % use CPU
    %num_segments = 1;
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
% Initialize everything in the specified segment, including both propagating forward and backward
% Everything is initialized because it's not sure which will be used.
% If forward propagation is done first, we need the backward data since in this code, both directions are, by default, considered.
% If backward propagation is done first, we need the forward data since in this code, both directions are, by default, considered.
[At_forward,Power_pump_forward,...
 N] = initialization(gain_rate_eqn,...
                     Nt,num_windows,...
                     num_zPoints,...
                     initial_condition);
T_delay = zeros(num_zPoints,1);

%% Propagations
% =========================================================================
% Start the first pulse propagation
% =========================================================================
gain_propagate();

%% Output:
if sim.gpu_yes
    [At_forward,Power_pump_forward,...
     N] = mygather2(1,num_zPoints,...
                    At_forward,Power_pump_forward,N);
end


% Transform them into arrays
At_out = cell(1,num_windows);
Power_out = cell(1,num_windows);
N_out = cell(1,num_windows);
for window_j = 1:num_windows
    At_out{window_j}   = struct('forward', cell2mat(At_forward (1,1,:,window_j)));

    Power_out{window_j} = struct('pump',struct('forward', cell2mat(Power_pump_forward (1,1,:,window_j))));

    N_out{window_j} = cell2mat(N(1,1,:,window_j));
end

%%
function gain_propagate()
%GAIN_PROPAGATE Runs the corresponding propagation method based on "direction".

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
        T_delay(1) = TCenter*dt(2);
    else
        T_delay(1) = 0;
    end
else
    T_delay(1) = 0;
end

% Put the first segment into GPU
if sim.gpu_yes
    [At_forward,Power_pump_forward,...
     N] = mygpuArray2(1,segments(1),...
                      At_forward,Power_pump_forward,...
                      N);
end

if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('SteppingCaller_rategain:ProgressBarNameError',...
              '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.0%%',sim.progress_bar_name),...
        'Name',sprintf('Running UPPE: %s...',sim.progress_bar_name),...
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
for ii = 2:num_zPoints
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('SteppingCaller_rategain:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end
    
    % =====================================================================
    % UPPE: Run the correct step function depending on the options chosen.
    % =====================================================================
    % Load the initial powers and field
    if ii == 2 % the first/starting z-index
        last_At_forward = At_forward(:,:,1,:); % = initial_condition.fields.forward
        last_Power_pump_forward = Power_pump_forward(:,:,1,:);
    end
    
    [last_At_forward,...
     last_Power_pump_forward,...
     N(:,:,ii-1,:)] = forward_stepping_RK4IP_rategain(last_At_forward,...
                                                      last_Power_pump_forward,...
                                                      N(:,:,ii-1,:),...
                                                      Nt,dt,num_windows,...
                                                      sim,gain_rate_eqn,...
                                                      SK_info,SRa_info,SRb_info,...
                                                      extended_haw,extended_hbw,...
                                                      extended_n2_prefactor,extended_expD,...
                                                      At_noise);
    
    % Update "forward" only
    At_forward(:,:,ii,:) = last_At_forward;
    Power_pump_forward(:,:,ii,:) = last_Power_pump_forward;
    
    % Save N
    if ii == num_zPoints % save the last N
        N(:,:,ii,:) = find_transient_N(sim,gain_rate_eqn,...
                                       Nt,dt,num_windows,...
                                       N(:,:,ii,:),...
                                       last_At_forward,...
                                       last_Power_pump_forward);

        if sim.gpu_yes
            N = mygather(N);
        end

        N = cellfun(@(x)x/gather(gain_rate_eqn.N_total),N,'UniformOutput',false);
    end
    
    % Whenever at the end of each segment, gather the current data
    % from GPU back to the RAM and then throw those in the next
    % segment to GPU.
    if rem(ii,zPoints_each_segment) == 0 || ii == num_zPoints
        ready_for_next_segment_into_GPU = true;
    else
        ready_for_next_segment_into_GPU = false;
    end
    
    % =====================================================================
    % Some post-stepping checks, saves, and updates
    % =====================================================================
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
            last_N_allWindows = cat(1,N{1,1,ii-1,:});

            for window_i = 1:num_windows
                last_Power_pump_forward{window_i} = interp1(t,last_Power_pump_forward_allWindows,t_each_window{window_i}+TCenter*dt(2),'linear','extrap');
                N{1,1,ii-1,window_i} = interp1(t,last_N_allWindows,t_each_window{window_i}+TCenter*dt(2),'linear','extrap');
                if sim.gpu_yes && mod(window_i,2) == 1
                    last_Power_pump_forward{window_i} = gather(last_Power_pump_forward{window_i});
                    N{1,1,ii-1,window_i} = gather(N{1,1,ii-1,window_i});
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
                T_delay(ii) = T_delay(ii-1) + gather(TCenter*dt(2));
            else
                T_delay(ii) = T_delay(ii-1) + TCenter*dt(2);
            end
        end
    end
    
    % Gather the current GPU data back to RAM and those in the next segment from RAM to GPU
    if ready_for_next_segment_into_GPU
        current_segment_idx = ceil(ii/zPoints_each_segment);
        next_segment_idx = current_segment_idx + 1;

        cumsum_segments = [0,cumsum(segments)];
        start_idx = cumsum_segments(current_segment_idx)+1;
        end_idx   = cumsum_segments(current_segment_idx+1);
        [At_forward,Power_pump_forward,...
         N] = mygather2(start_idx,end_idx,...
                        At_forward,Power_pump_forward,N);

        if next_segment_idx > 0 && next_segment_idx <= length(segments)
            start_idx = cumsum_segments(next_segment_idx)+1;
            end_idx   = cumsum_segments(next_segment_idx+1);
            [At_forward,Power_pump_forward,...
             N] = mygpuArray2(start_idx,end_idx,...
                              At_forward,Power_pump_forward,N);
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

end

end

%% initialization
function [At_forward,Power_pump_forward,...
          N] = initialization(gain_rate_eqn,...
                              Nt,num_windows,...
                              zPoints,...
                              initial_condition)
%INITIALIZATION initializes "signal_fields" and "Powers"
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
% Pump and ASE
Power_pump_forward = initialize_zeros();
At_forward = initialize_zeros();
N = initialize_zeros();

% -------------------------------------------------------------------------
% Put in the necessary information
% -------------------------------------------------------------------------
% Pump power
for window_i = 1:num_windows
    At_forward{1,1,1,window_i} = initial_condition(window_i).fields.forward;
    
    Power_pump_forward{1,1,1,window_i} = gain_rate_eqn.copump_power*ones(Nt,1);

    num_levels = size(initial_condition(window_i).population,2); % this is actually num_levels-1
    if size(initial_condition(window_i).population,1) == 1
        N(1,1,:,window_i) = mat2cell(repmat(initial_condition(window_i).population*gain_rate_eqn.N_total,Nt,1,1),Nt,num_levels,ones(1,zPoints)); % size: (1,num_levels-1)
    else
        N(1,1,:,window_i) = mat2cell(       initial_condition(window_i).population*gain_rate_eqn.N_total,        Nt,num_levels,ones(1,zPoints)); % size: (1,num_levels-1)
    end
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
    varargout{i}(:,:,starti:endi,:) = cellfun(@gather,varargin{i}(:,:,starti:endi,:),'UniformOutput',false);
end

end

%% CLEANMEUP
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end