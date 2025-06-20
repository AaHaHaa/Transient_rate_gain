function foutput = UPPE_propagate_without_adaptive(fiber,initial_condition,sim,gain_rate_eqn)
%UPPE_PROPAGATE_WITHOUT_ADAPTIVE Propagate an initial pulse 
%through an arbitrary distance of an optical fiber without the use of 
%adaptive-step method.
%
% -------------------------------------------------------------------------
%
%   "fiber" is a structure with the fields:
%
%       Basic properties -->
%
%           betas - a (?,1) array
%                   betas(i) = (i-1)th order dispersion coefficient for each mode, in ps^n/m
%
%           n2 - the nonlinear coefficient (default to 2.3e-20 if not set)
%
%           SR - SR tensor, in m^-2
%           L0 - length of fiber, in m
%           material - 'silica', 'chalcogenide', or 'ZBLAN' (default: 'silica')
%                      This is used to determine the Raman parameters to use.
%
% -------------------------------------------------------------------------
%
%   "initial_condition" is a 1-by-2 structure array, each with the fields:
%
%       dt - time step, in ps
%       fields.forward - initial forward field in time domain, in W^1/2; (N,1).
%       fields.backward - initial backward field in time domain, in W^1/2; (N,1).
%
% -------------------------------------------------------------------------
%
%   "sim" is a structure with the fields:
%
%       Basic settings -->
%
%           betas - the betas for the slowly varying approximation and the moving frame, 
%                   that is to say, fiber.betas([1 2],:) = fiber.betas([1 2],:) - sim.betas;
%                   (2,1) column vector;
%                   if not set, no "sim.betas", the simulation will be run relative to the first mode
%           f0 - center frequency, in THz
%           dz - step size, in m
%           save_period - spatial period between saves, in m
%                         0 = only save input and output (save_period = fiber.L0)
%
%       Algorithms to use -->
%
%           gpu_yes - 1(true) = GPU
%                     0(false) = CPU
%
%                     Whether or not to use the GPU. Using the GPU is HIGHLY recommended, as a speedup of 50-100x should be possible.
%
%           include_Raman - 0(false) = ignore Raman effect
%                           1(true) = Either (a) Raman model approximated analytically by a single vibrational frequency of silica molecules
%                                                (Ch. 2.3, p.42, Nonlinear Fiber Optics (5th), Agrawal)
%                                                Extensions to other materials are also included. Please check Raman_model().
%                                     or     (b) Raman model including the anisotropic contribution
%                                                ("Ch. 2.3, p.43" and "Ch. 8.5, p.340", Nonlinear Fiber Optics (5th), Agrawal)
%                                                For more details, please read "Raman response function for silica fibers", by Q. Lin and Govind P. Agrawal (2006)
%
%       Scaled Fourier transform -->
%
%           cs - scaled ratio (an integer) for narrowband transformation.
%                It must be a divisor of Nt.
%
%       Others -->
%
%           pulse_centering - 1(true) = center the pulse according to the time window, 0(false) = do not
%                             The time delay will be stored in time_delay after running Transient_gain_UPPE_propagate().
%           gpuDevice.Index - a scalar; the GPU to use
%           gpuDevice.Device - the output of MATLAB "gpuDevice(gpu_index)"
%           cuda_dir_path - path to the cuda directory into which ptx files will be compiled and stored
%           progress_bar - 1(true) = show progress bar, 0(false) = do not
%                          It'll slow down the code slightly. Turn it off for performance.
%           progress_bar_name - the name of the UPPE propagation shown on the progress bar.
%                               If not set (no "sim.progress_bar_name"), it uses a default empty string, ''.
%
% Output:
%   foutput is a (1,num_windows) structure array with fields:
%
%       z - the propagation length of the saved points
%       dz - the step size at each saved points
%       fields.forward - (Nt, 1, num_save_points) array with the forward field at each save point
%       fields.backward - (Nt, 1, num_save_points) array with the backward field at each save point
%       dt - time grid point spacing, to fully identify the field
%       betas - the [betas0;betas1] used for the moving frame
%       seconds - time spent in the main loop
%       t_delay - the time delay of each pulse which is centered in the time window during propagation
%       Power.pump.forward  - (Nt,1,num_save_points); the forward pump power along the fiber
%       Power.pump.backward - (Nt,1,num_save_points); the backward pump power along the fiber
%       population - (Nt,num_levels,num_save_points); the doped ion density of upper levels
%       shot_noise - a (1,num_windows) cell array, each with a (Nt,1) array of shot noise in its window

%% Check the validity of input parameters
if sim.gpu_yes
    try
        gpuDevice;
    catch
        error('UPPE_propagate_without_adaptive:GPUError',...
              'No GPU is detected. Please set "sim.gpu_yes=false".');
    end
end

num_windows = length(initial_condition);
if mod(num_windows,2) ~= 0
    error('UPPE_propagate_with_adaptive:num_windowsError',...
          'There must be even number of time windows, alternating between with and without coherent fields.');
end

Nt = arrayfun(@(x)size(x.forward,1),[initial_condition.fields]);
dt = [initial_condition.dt];
if any(Nt-mean(Nt)) || any(abs(dt(2:2:end)-mean(dt(2:2:end))) > eps(mean(dt(2:2:end)))*3) % dt is double-precison whose computation can have errors, so we need eps() to check its equality
    error('UPPE_propagate_with_adaptive:NtError',...
          ['(1) All Nt need to be the same but dt can be different.\n',...
           '(2) All Nt and dt for coherent fields need to be the same.']);
end
Nt = Nt(1); % because Nt is the same for all windows, we take only the first one

for window_i = 1:num_windows
    if size(initial_condition(window_i).fields.forward,3) > 1 || size(initial_condition(window_i).fields.backward,3) > 1
        error('UPPE_propagate_with_adaptive:initial_conditionError',...
              'Make sure that there is only one field in each direction and temporal window of "initial_condition".');
    end
end

if sim.save_period == 0
    sim.save_period = fiber.L0;
end

% MPA performs M "sim.dz/M" in parallel
num_zSteps = fiber.L0/sim.dz;
num_zPoints_persave = sim.save_period/sim.dz;
num_saveSteps = fiber.L0/sim.save_period;

% Because of machine error, "rem" alone to determine divisibility isn't
% enough. eps() is included here.
if rem(num_zSteps,1) && rem(num_zSteps+eps(num_zSteps),1) && rem(num_zSteps-eps(num_zSteps),1)
    error('UPPE_propagate:SizeIncommemsurateError',...
        'The step size is %f m and the fiber length is %f m, which are not commensurate', sim.dz, fiber.L0)
else
    num_zSteps = round(num_zSteps);
end

if rem(num_zPoints_persave,1) && rem(num_zPoints_persave+eps(num_zPoints_persave),1) && rem(num_zPoints_persave-eps(num_zPoints_persave),1)
    error('UPPE_propagate:SizeIncommemsurateError',...
        'The step size is %f m and the save period is %f m, which are not commensurate', sim.dz, sim.save_period)
else
    num_zPoints_persave = round(num_zPoints_persave);
end

if rem(num_saveSteps,1) && rem(num_saveSteps+eps(num_saveSteps),1) && rem(num_saveSteps-eps(num_saveSteps),1)
    error('UPPE_propagate:SizeIncommemsurateError',...
        'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saveSteps = round(num_saveSteps);
end

%% the scaled factor of the scaled Fourier transform
cs = dt/sim.incoherent.dt;
cs(2:2:end) = sim.cs;

%% Pre-calculate the dispersion term
% The "Omega" here is offset frequency (omega - omega0), omega: true angular frequency
%                                                        omega0: central angular frequency (=2*pi*f0)
Omega = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/(Nt*dt(2)); % in rad/ps, in the order that the ifft gives

% The dispersion term in the UPPE, in frequency space
peak_power = zeros(1,num_windows);
for window_i = 1:num_windows
    peak_power(window_i) = max(abs(initial_condition(window_i).fields.forward).^2,[],1);
end
[~,max_idx] = max(peak_power);
dominant_field = initial_condition(max_idx).fields.forward;
[D_op,sim] = calc_D_op(fiber,sim,Nt,dt(2),Omega,dominant_field);

%% Create a damped frequency window to prevent aliasing
sim.damped_freq_window = create_damped_freq_window(Nt);

%% Pre-calculate the factor used in UPPE (the nonlinear constant)
c = 2.99792458e-4; % speed of ligth m/ps
if ~isfield(fiber,'n2') || isempty(fiber.n2)
    fiber.n2 = 2.3e-20; % m^2/W
end
n2_prefactor = 1i*fiber.n2*(Omega+2*pi*sim.f0)/c; % m/W

% Incorporate the damped window to the n2 term to remove generation of
% frequency component near the window edge.
n2_prefactor = n2_prefactor.*sim.damped_freq_window;

% Scaled nonlinearity due to the narrowband transformation (scaled Fourier transform)
if sim.cs > 1
    n2_prefactor = n2_prefactor/sim.cs;
end

%% Extend the time window for the acyclic-convolution operation
acyclic_conv_ratio = 2;
gain_rate_eqn.acyclic_conv_stretch = @(x) acyclic_conv_ratio*x-1;
gain_rate_eqn.acyclic_conv_shrink = @(y) (y+1)/acyclic_conv_ratio; % reverse of acyclic_conv_stretch(); this is used to remove population in the extended time window in solve_Power() in solve_gain_rate_eqn_helpers()

original_f = linspace(-floor(Nt/2), floor((Nt-1)/2), Nt)'/Nt;
extended_f = linspace(-floor(gain_rate_eqn.acyclic_conv_stretch(Nt)/2), floor((gain_rate_eqn.acyclic_conv_stretch(Nt)-1)/2), gain_rate_eqn.acyclic_conv_stretch(Nt))'/gain_rate_eqn.acyclic_conv_stretch(Nt);

extended_D_op = ifftshift(interp1(original_f,fftshift(D_op,1),extended_f,'linear','extrap'),1);
extended_n2_prefactor = ifftshift(interp1(original_f,fftshift(n2_prefactor,1),extended_f,'linear','extrap'),1);
for window_i = 1:2
    gain_rate_eqn.extended_cross_sections{window_i} = ifftshift(interp1(original_f,fftshift(gain_rate_eqn.cross_sections{window_i},1),extended_f,'linear','extrap'),1);
    gain_rate_eqn.extended_E_photon{window_i} = ifftshift(interp1(original_f,fftshift(gain_rate_eqn.E_photon{window_i},1),extended_f,'linear','extrap'),1); % J
    % Just a check below:
    gain_rate_eqn.extended_cross_sections{window_i}(gain_rate_eqn.extended_cross_sections{window_i}<0) = 0;
    gain_rate_eqn.extended_E_photon{window_i}(gain_rate_eqn.extended_E_photon{window_i}<0) = 0;
end

gain_rate_eqn.acyclic_zero_padding = @(x) [x;zeros(gain_rate_eqn.acyclic_conv_stretch(length(x))-length(x),size(x,2))];

%% Pre-compute the dispersion operator
extended_expD = exp(extended_D_op*sim.dz/2);

%% Set up the GPU details
if sim.gpu_yes
    [sim.gpuDevice.Device,...
     sim.cuda_SRSK,sim.cuda_num_operations_SRSK] = setup_stepping_kernel(sim,gain_rate_eqn.acyclic_conv_stretch(Nt));
end

%% Pre-compute the Raman response in frequency space
if ~isfield(fiber,'material')
    fiber.material = 'silica';
end
[fiber,extended_haw,extended_hbw] = Raman_model(fiber,sim,gain_rate_eqn.acyclic_conv_stretch(Nt),dt(2));

%% Work out the overlap tensor details
[SK_info,SRa_info,SRb_info] = calc_SRSK(fiber,sim,1);

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saveSteps + 1);
num_zPoints = num_zSteps + 1;

%% Put data in GPU if using GPU
if sim.gpu_yes
    for window_i = 1:2
        gain_rate_eqn.extended_cross_sections{window_i} = gpuArray(gain_rate_eqn.extended_cross_sections{window_i});
        gain_rate_eqn.extended_E_photon{window_i} = gpuArray(gain_rate_eqn.extended_E_photon{window_i});
    end
    
    extended_n2_prefactor = gpuArray(extended_n2_prefactor);
    extended_expD = gpuArray(extended_expD);

    extended_haw = gpuArray(extended_haw);
    extended_hbw = gpuArray(extended_hbw);
    
    for window_i = 1:num_windows
        initial_condition(window_i).fields.forward = gpuArray(initial_condition(window_i).fields.forward);
        initial_condition(window_i).fields.backward = gpuArray(initial_condition(window_i).fields.backward);
    end
end

%% Modified shot noise for noise modeling
At_noise = cell(2,num_windows); % 2 for forward and backward
for direction_i = 1:2
    for window_i = 1:num_windows
        if mod(window_i,2) == 1
            At_noise{direction_i,window_i} = shot_noise(gain_rate_eqn.acyclic_conv_stretch(Nt),dt(window_i),sim.incoherent.f0,cs(window_i));
        else
            At_noise{direction_i,window_i} = shot_noise(gain_rate_eqn.acyclic_conv_stretch(Nt),dt(window_i),sim.f0,cs(window_i));
        end

        if sim.gpu_yes
            At_noise{direction_i,window_i} = gpuArray(At_noise{direction_i,window_i});
        end
    end
end

%% Run the step function over each step
run_start = tic;
% -------------------------------------------------------------------------
[A_out,Power,...
 T_delay_out,...
 N] = SteppingCaller_rategain(sim,gain_rate_eqn,...
                              num_zPoints,save_points,num_zPoints_persave,...
                              initial_condition,...
                              extended_n2_prefactor,...
                              SK_info,SRa_info,SRb_info,...
                              extended_expD,...
                              extended_haw,extended_hbw,...
                              At_noise);

% -------------------------------------------------------------------------
% Just to get an accurate timing, wait before recording the time
if sim.gpu_yes
    sim.betas = gather(sim.betas);
    for window_i = 1:num_windows
        At_noise{1,window_i} = gather(At_noise{1,window_i});
    end
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Shot noise
At_noise = At_noise(1,:);
for window_i = 1:num_windows
    Aw_noise = ifft(At_noise{window_i});
    amp_Aw_noise = abs(Aw_noise);
    ang_Aw_noise = unwrap(angle(Aw_noise));
    amp_Aw_noise = ifftshift(interp1(extended_f,fftshift(amp_Aw_noise,1),original_f,'linear','extrap'),1);
    ang_Aw_noise = ifftshift(interp1(extended_f,fftshift(ang_Aw_noise,1),original_f,'linear','extrap'),1);
    Aw_noise = amp_Aw_noise.*exp(1i*ang_Aw_noise);
    At_noise{window_i} = fft(Aw_noise,[],1);
end

%% Save the results in a struct
foutput = struct('z', sim.save_period*(0:num_saveSteps)',...
                 'fields', A_out,...
                 'dt', num2cell(dt),...
                 'betas', sim.betas,...
                 'seconds', fulltime,...
                 't_delay', T_delay_out,...
                 'Power', Power,...
                 'population', N,...
                 'shot_noise', At_noise);

end