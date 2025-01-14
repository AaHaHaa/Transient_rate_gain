function [fiber,sim] = load_default_UPPE_propagate( input_fiber,input_sim )
%LOAD_DEFAULT_UPPE_PROPAGATE It loads the default settings for "fiber"
%and "sim" for different types of modes used.
%
%   If a user has specified some of the parameters of "fiber" and "sim",
%   user-defined one will be chosen instead of the default ones.
%
%   If you want to use a different one, specify it as fiber.xxx or sim.xxx 
%   and send them into this function. The output will use your parameters 
%   besides other default parameters a user doesn't specify.
%
%% ========================================================================
%   Because some parameters are correlated to each other, sometimes it's 
%   necessary to load input parameters first.
%   However, this makes loading default parameters a complicated function.
%
%   The procedure of loading default parameters is described below.
%   The order matters!
% -------------------------------------------------------------------------
%
%   <-- Uncorrelated parameters are loaded directly -->
%
%       sim.f0 - depend on input f0 or lambda0
%                If no input f0 or lambda0, f0=3e5/1030e-9 (THz)
%
% If there's user-defined one, use user's for the parameters below.
% Below I list the default values -- >
%
%       fiber.material = 'silica';
%       fiber.n2 = 2.3e-20;
%
%       sim.dz = 1000e-6;
%       sim.save_period = 0;
%
%       sim.scalar = true;
%
%       sim.adaptive_dz.threshold = 1e-6;
%
%       sim.gpu_yes = true;
%       sim.include_Raman = true;
%
%       sim.pulse_centering = true;
%       sim.gpuDevice.Index = 1;
%       sim.progress_bar = true;
%       sim.progress_bar_name = '';
%       sim.cuda_dir_path = [folder_of_this_function 'UPPE/cuda'];
%
%   <-- Correlated parameters are loaded based on the input or default -->
%
%   single-mode -- >
%
%       sim.midx = [] (not used)
%
%       Assume 1030 nm for positive dispersion if lambda0 < 1300 nm (~ZDW for a silica fiber),
%           fiber.betas = [8.867e6; 4.903e3; 0.0208; 33.3e-6; -27.7e-9];
%           fiber.MFD = 6.2; % um; 1030nm from Thorlabs 1060XP
%       Assume 1550 nm for negative dispersion if lambda0 > 1300 nm (~ZDW for a silica fiber),
%           fiber.betas = [5.87e6; 4.91e3; -0.0167; 1.04e-6; -324e-9];
%           fiber.MFD = 9.5; % um; 1030nm from Thorlabs 1060XP
%
%       fiber.SR = (1) input SR, if there's input SR (input SR precedes over input MFD)
%                  (2) 1/Aeff, if (a) there's input MFD
%                                 (b) MFD is taken from the default one and there's no input MFD
%
%   fiber.L0 = (1) input L0
%              (2) 2 (m)
%
%% ========================================================================
%   Details of each parameter
% -------------------------------------------------------------------------
%
% Example Use:
%
%    % User-defined parameters
%    fiber.betas = [0 0 0.02 0];
%    fiber.L0 = 3;
%    
%    % Incorporate default settings
%    [fiber,sim] = load_default_UPPE_propagate(fiber,[]); % single_mode
%
%    % If there are "sim" settings
%    sim.adaptive_dz.model = 0;
%    [fiber,sim] =  load_default_UPPE_propagate(fiber,sim); % single_mode
%
%    % Use only user-defined "sim", not "fiber"
%    [fiber,sim] = load_default_UPPE_propagate([],sim); % single_mode
%
% -------------------------------------------------------------------------
%
%	Additional parameters:
%
%       input_sim.lambda0 - central wavelength, in m
%       input_fiber.MFD - mode field diameter for calculating Aeff for SR in single mode; only used in single mode, in um
%
% -------------------------------------------------------------------------
%
%   If both "lambda0" and "f0" are set by users, the final value will depend on "f0".
%
% -------------------------------------------------------------------------
%
%   "fiber" is a structure with the fields:
%
%       Basic properties -->
%
%           betas - a (?,nm) matrix; "nm" = num_spatial_modes if under scalar fields;
%                                    otherwise, "nm" can be both num_spatial_modes or 2*num_spatial_modes depending on whether there's birefringence.
%                   betas(i, :) = (i-1)th order dispersion coefficient for each mode, in ps^n/m
%
%           n2 - the nonlinear coefficient (default to 2.3e-20 if not set)
%
%           SR - SR tensor, in m^-2
%           L0 - length of fiber, in m
%           fiber_type - 'silica', 'chalcogenide', or 'ZBLAN' (default: 'silica')
%
% -------------------------------------------------------------------------
%
%   "initial_condition" is a structure with the fields:
%
%       dt - time step
%       fields.forward - initial forward field, in W^1/2, (N,1).
%       fields.backward - initial backward field, in W^1/2, (N,1).
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
%       Adaptive method -->
%
%           adaptive_dz.threshold - a scalar;
%                                       the accuracy used to determined whether to increase or decrease the step size.
%           adaptive_dz.max_dz - a scalar; the maximum adaptive step size
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
%       Others -->
%
%           pulse_centering - 1(true) = center the pulse according to the time window, 0(false) = do not
%                             The time delay will be stored in time_delay after running UPPE_propagate().
%           gpuDevice.Index - a scalar; the GPU to use
%           gpuDevice.Device - the output of MATLAB "gpuDevice(gpu_index)"
%           cuda_dir_path - path to the cuda directory into which ptx files will be compiled and stored
%           progress_bar - 1(true) = show progress bar, 0(false) = do not
%                          It'll slow down the code slightly. Turn it off for performance.
%           progress_bar_name - the name of the UPPE propagation shown on the progress bar.
%                               If not set (no "sim.progress_bar_name"), it uses a default empty string, ''.
%
% =========================================================================

%% Current path (or the folder where this "load_default_UPPE_propagate.m" is)
if ispc
    sep = '\';
else % unix
    sep = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep);
upper_folder = current_path(1:sep_pos(end-1));

%% Default settings below:

% Supperss warnings generated by the function, "catstruct", due to there
% are definitely duplicate elements between default and input.
warning('off','catstruct:DuplicatesFound');

if ~exist('input_fiber','var')
    input_fiber = [];
end
if ~exist('input_sim','var')
    input_sim = [];
end

% -------------------------------------------------------------------------
% Set some default parameters early here because these parameters will be
% used for loading files if multimode.
% If single-mode, only "fiber.lambda0" or "fiber.f0" is important.
%
% Set lambda0 below:
% Priority: user_f0 > user_lambda0 > default_lambda0
% This lambda0 will then be set to f0 later for the UPPE functions to use.
% -------------------------------------------------------------------------
c = 2.99792458e-4; % speed of ligth, m/ps

% Get lambda0 from input f0 or lambda0
if isfield(input_sim,'f0')
    default_sim.lambda0 = c/input_sim.f0;
else
    if isfield(input_sim,'lambda0')
        default_sim.lambda0 = input_sim.lambda0;
    else
        default_sim.lambda0 = 1030e-9;
    end
end

% -------------------------------------------------------------------------
% fiber
% -------------------------------------------------------------------------
% Basic properties
% (Load betas, SR tensors, and Aeff below)
default_fiber.material = 'silica';
default_fiber.n2 = 2.3e-20; % m^2 W^-1
default_sim.midx = 1;

if default_sim.lambda0 > 1e-6 && default_sim.lambda0 < 1.3e-6 % lambda(zero dispersion)=1.3um; assume 1030nm for normal dispersion
    default_fiber.betas = [8.8268e6; 4.8821e3; 0.0209; 32.9e-6; -26.7e-9];
    default_fiber.MFD = 5.95; % um; 1030nm from Thorlabs's 1060XP
elseif default_sim.lambda0 < 1e-6 % assume at 920 nm, where Nd is used
    default_fiber.betas = [9.8810e6; 4.8872e3; 0.0315; 2.0457e-5; 1.2737e-9];
    default_fiber.MFD = 4; % um; 1030nm from IXblue's IXF-2CF-PAS-PM-4-80-0.16-P
else % assume 1550nm for anomalous dispersion
    default_fiber.betas = [5.8339e6; 4.8775e3; -0.0123; 0.1049e-6; -378.3e-9];
    default_fiber.MFD = 8.09; % um; 1550nm from Thorlabs's 1060XP
end
% Load "input_fiber" into "default_fiber" first to calculate SR.
if isfield(input_fiber,'MFD')
    default_fiber.MFD = input_fiber.MFD;
end

Aeff = pi*(default_fiber.MFD/2)^2*1e-12; % m^2; effective area of the SMF
default_fiber.SR = 1/Aeff;

% Because MFD is used to calculate SR for single-mode case,
% clear it if input_fiber has SR already.
if isfield(input_fiber,'SR')
    Aeff = 1./input_fiber.SR;
    default_fiber.MFD = sqrt(Aeff/pi)*2;
end

% -------------------------------------------------------------------------
% sim
% -------------------------------------------------------------------------
% Basic settings
default_sim.f0 = c/default_sim.lambda0; % THz
default_sim.dz = 1000e-6; % m
default_sim.save_period = 0; % m

% Adaptive method
% Threshold error for adaptive RK4IP
default_sim.adaptive_dz.threshold = 1e-7; % the threshold of the adaptive method
                                              % Recommended value is less than 1e-5.
                                              % Values larger than 1e-3 are too large.

% Algorithms to use
default_sim.gpu_yes = true;
default_sim.include_Raman = true; % consider Raman

% Others
default_sim.pulse_centering = true; % center the pulse according to the time window
default_sim.gpuDevice.Index = 1; % the gpuDevice to use
default_sim.progress_bar = true;
default_sim.progress_bar_name = '';
default_sim.cuda_dir_path = fullfile(upper_folder,'cuda');

%%
% =========================================================================
% Merge settings with the input, which have higher priorities than the
% default ones.
% =========================================================================
if isempty(input_fiber)
    fiber = default_fiber;
elseif isstruct(input_fiber)
    fiber = catstruct(default_fiber, input_fiber);
else
    error('LoadDefaultUPPEPropagate:InputFiberError',...
            '"input_fiber" should be a "structure".');
end
if isempty(input_sim)
    sim = default_sim;
elseif isstruct(input_sim)
    sim = catstruct(default_sim, input_sim);
else
    error('LoadDefaultUPPEPropagate:InputSimError',...
            '"input_sim" should be a "structure".');
end

end