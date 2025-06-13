function gain_rate_eqn = gain_info( sim,gain_rate_eqn,lambda )
%GAIN_INFO Computes several information related to the gain medium.
%
% =========================================================================
% =============== Call this function with the following code ==============
% =========================================================================
% f = ifftshift( (-N/2:N/2-1)/N/dt + sim.f0 ); % in the order of "Omega" in the Transient_gain_UPPE_propagate()
% c = 299792.458; % nm/ps
% lambda = c./f; % nm
%
% gain_rate_eqn = gain_info( sim,gain_rate_eqn,lambda );
% output_field = Transient_gain_UPPE_propagate(fiber,input_field,sim,gain_rate_eqn);
%
% =========================================================================
%
%   gain_rate_eqn:
%
%       fiber info -->
%
%           core_diameter - um; where the doped ion stays
%           cladding_diameter - um
%           core_NA - numerical aperture of the gain fiber
%
%       doped ion info -->
%
%           absorption_wavelength_to_get_N_total - nm
%           absorption_to_get_N_total - dB/m
%
%       pump info -->
%
%           pump_wavelength - nm
%           pump_direction - "Don't set this parameter" since it'll be automatically determined from "copump_power" and "counterpump_power" below.
%                            'co', 'counter', or 'bi'
%                            'co' reads only the copump_power, 'counter' reads only the counterpump_power, while 'bi' reads both.
%           copump_power - W
%           counterpump_power - W
%
%       ASE info -->
%
%           sponASE_spatial_modes - the number of ASE's supported spatial modes
%                                   In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE.
%                                   If empty like [], it's length(sim.midx).
%
%       rate equation model algorithm info -->
%
%           ignore_ASE - 1(true) or 0(false)
%           max_iterations - the maximum number of iterations
%           tol - the tolerance of this iteration loop. If the difference between the last two results is smaller than this tolerance, we're done. 
%           verbose - show the information(final pulse energy) during iterations
%           memory_limit: the RAM limit for this computation
%                         (This can be found by default.)
%
%   lambda - the wavelengths of the computational region; in "nm"
%
%   =======================================================================
%   Output (in "gain_rate_eqn"):
%
%       cross_sections_pump - um^2
%       cross_sections - um^2
%       overlap_factor - 1/um^2
%       N_total - a scalar; in "1/um^3"
%       ASE_factor - num_polarizations*sponASE_spatial_modes
%       E_photon - h*nu of each window; in "J"

%% Add the folder of functions of gain-rate-equation model and its functions
% Besides loading mode coupling folder, this "sep_char" is also used in GPU setup below.
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));
addpath([upper_folder 'Gain_rate_eqn/'],...
        [upper_folder 'Gain_rate_eqn/gain cross sections/'],...
        [upper_folder 'Gain_rate_eqn/Judd-Ofelt theory/Material data/']);

%% Add more-verbose parameters
% Activate ASE computation under counterpumping or bidirectional pumping.
% Since back-and-forth iterations are needed for counterpumping or bidirectional pumping anyway, it'll just run with ASE.
if gain_rate_eqn.counterpump_power ~= 0
    gain_rate_eqn.ignore_ASE = false;
end

gain_rate_eqn.include_ASE = ~gain_rate_eqn.ignore_ASE;

%% "lambda" error check
% "lambda" must be a column vector.
for window_i = 1:2
    if size(lambda{window_i},1) == 1
        lambda{window_i} = lambda{window_i}.';
    end
    if any(lambda{window_i}(:)<0)
        error('gain_info:lambdaError',...
              ['Wavelength, the "lambda" input variable, cannot be negative.\n',...
               'If possible, don''t use the pulse frequency as the center of the frequency window, which helps offset the frequency window from getting close to zero.\n',...
               'Use find_tw_f0() to help do this.']);
    end
    % MATLAB before 2017 doesn't have "issortedrows()" and 'monotonic' argument in "issorted()"
    MATLAB_version = version('-release'); MATLAB_version = str2double(MATLAB_version(1:4));
    if MATLAB_version < 2017
        do_ifftshift_on_lambda = (issorted(lambda{window_i}) || issorted(flipud(lambda{window_i})));
    else
        do_ifftshift_on_lambda = issortedrows(lambda{window_i},'monotonic'); % = issorted(lambda,'monotonic');
    end
    if do_ifftshift_on_lambda
        lambda{window_i} = ifftshift(lambda{window_i},1);
    end
end

%% Narrowband transformation of the coherent fields due to the scaled Fourier transform
% lambda needs to be re-defined
if sim.cs > 1
    for window_i = 1:2
        lambda{window_i} = lambda{window_i}(1:sim.cs:end);
    end
end

%% Function container to read necessary parameters based on the gain medium to use later
gain_func = gain_medium();

% Read necessary parameters and cross sections based on the gain medium to use
gain_rate_eqn = gain_func.load_medium_parameters(gain_rate_eqn);
cross_sections = cell(1,2);
for window_i = 1:2
    [gain_rate_eqn,...
     cross_sections{window_i},cross_sections_pump,...
     GSA_find_Ntotal] = gain_func.load_cross_sections(gain_rate_eqn,lambda{window_i});
end

%% Overlap factor of the field and the doping area
% pump overlap factor = A_doping/A_cladding
overlap_factor.pump = 1/(pi*(gain_rate_eqn.cladding_diameter/2)^2);

% signal overlap factor
V = 2*pi*(gain_rate_eqn.core_diameter/2)/(lambda{2}(1)*1e-3)*gain_rate_eqn.core_NA; % V-number, normalized frequency
if V < 0.8 || V > 2.8
    warning('For the computation of the fundamental mode, I use "Whitley''s Gaussian-mode approximation". It works only in the range of V=0.8-2.8.\nCurrent V number is %4.2f',V);
end
% w_over_a = 0.65+1.619*V^(-3/2)+2.879*V^(-6); % Marcuse et al., "Loss analysis of single-mode fiber splices" (1977)
w_over_a = 0.616+1.66*V^(-3/2)+0.987*V^(-6); % Whitley et al., "Alternative Gaussian spot size polynomial for use with doped fiber amplifiers" (1993)
overlap_factor.signal = ( 1-exp(-2/w_over_a^2) )/(pi*(gain_rate_eqn.core_diameter/2)^2); % 1/um^2

%% Doped ion density
% For small-signal absorption, N1~0 and N0~N_total.
% pump is proportional to "exp(-integral2(overlap_factor.pump)*N_total*absorption_cross_section*L)", L: propagation length
% absorption_dB/m =10*log10( exp(-integral2(overlap_factor.pump)*N_total*absorption_cross_section*(L=1m)) )
N_total = log(10^(gain_rate_eqn.absorption_to_get_N_total/10))./(((gain_rate_eqn.core_diameter/gain_rate_eqn.cladding_diameter)^2)*GSA_find_Ntotal*1e6); % doped ion density based on absorption at a specific wavelength; in "1/um^3"

%% Read necessary parameters based on the gain medium to use
gain_rate_eqn = gain_func.load_N_related_parameters(gain_rate_eqn,N_total);

%% Check the validity of the code
% Number of ASE spatial modes
% For LMA fibers, there are more than one spatial modes for ASE although
% the signal field mostly stays only within the fundamental mode. In this
% situation, the simulations are mostly run with single mode, so
% consideration of multimode ASE is included with this 
% sponASE_spatial_modes factor.
if gain_rate_eqn.include_ASE
    if isempty(gain_rate_eqn.sponASE_spatial_modes)
        gain_rate_eqn.sponASE_spatial_modes = 1; % only single-mode is supported now
    else
        if gain_rate_eqn.sponASE_spatial_modes < 1 % only single-mode is supported now
            error('gain_info:NumASEError',...
                  'The number of ASE spatial modes need to be larger than that of the signal field.');
        end
    end
else
    gain_rate_eqn.sponASE_spatial_modes = 0;
end

%% Query the gpuDevice
if sim.gpu_yes
    if ~isfield(sim,'gpuDevice') || ~isfield(sim.gpuDevice,'Device')
        sim.gpuDevice.Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
    end
end

%% Find the memory limit
if ~isfield(gain_rate_eqn,'memory_limit')
    if sim.gpu_yes
        gain_rate_eqn.memory_limit = sim.gpuDevice.Device.AvailableMemory/2;
    else
        if ispc % Windows
            userview = memory;
            gain_rate_eqn.memory_limit = userview.MaxPossibleArrayBytes/2; % B
        elseif isunix % Unix, linux
            [~,w] = unix('free -b | grep Mem'); % Display the memory info in Bytes
            stats = str2double(regexp(w, '[0-9]*', 'match'));
            %memsize = stats(1)/1e6;
            freemem = stats(end); % B; availabel memory
            gain_rate_eqn.memory_limit = freemem/2;
        else % iOS
            error('gain_info:OSError',...
                  'iOS is not well supported yet.');
        end
    end
end

%% Reorganize cross sections into an array for faster numerical computations
% Dimension: (Nt,num_cross_sections)
cross_sections_pump = cell2mat(struct2cell(cross_sections_pump)).'; % size: (1,num_cross_sections)

for window_i = 1:2
    cross_sections{window_i} = structfun(@(x) x.',cross_sections{window_i},'UniformOutput',false); % change it to the size (1,Nt) for each cross-section data
    cross_sections{window_i} = cell2mat(struct2cell(cross_sections{window_i})).'; % size: (Nt,num_cross_sections)
end

%% Population parameters
gain_rate_eqn.N.MPA = struct('tol', 1e-20,...
                             'min_iter', 2,...
                             'max_iter', 200);

% Maximum number of iterations for MPA
% Since I realized that MPA can always converge, even with a super long
% time window, for a two-level system (such as Yb), I set the maximum
% number of iterations to start considering switching to the 
% time-sequential method to be a huge value simply to let MPA run until it 
% converges. 
% For example, for a 100-Hz Yb-doped fiber amplifier (without ASE effect), 
% MPA needs up to ~110 iterations to converge.
% In principle, even up to this massive iterations, MPA should still be
% faster than slow sequential temporal computations with a for-loop.
% However, I realized that for a multi-level system, population, if not
% converged within a few iterations, will simply diverge. Therefore, its
% maximum number of iterations to start considering switching to the 
% time-sequential method is set smaller to trigger the sequential 
% computation ASAP.
if length(gain_rate_eqn.energy_levels) == 2 % two-level system (such as Yb)
    gain_rate_eqn.N.MPA.start_sequential_factor = 2;
else % multi-level system
    gain_rate_eqn.N.MPA.start_sequential_factor = 20;
end

% Below sets the number meaning that the tolerance must be reduced by 
% "tol_ratio" times as the initial tolerance during finding N(t) under the 
% periodic boundary condition.
% Check find_periodic_transient_N() for details.
gain_rate_eqn.N.periodic.tolf = 0.01;

%% ASE factor
% Even in scalar computations, the gain is saturated by the spontaneous emission of two polarizations
num_polarizations = 2;
ASE_factor = num_polarizations*gain_rate_eqn.sponASE_spatial_modes;

%% Photon energy
h = 6.62607015e-34;
c = 299792.458e12;
E_photon = cell(1,2);
for window_i = 1:2
    E_photon{window_i} = h*(c./lambda{window_i});
end

%% Put all information into "gain_rate_eqn"
gain_rate_eqn.cross_sections_pump = cross_sections_pump;
gain_rate_eqn.cross_sections = cross_sections;
gain_rate_eqn.overlap_factor = overlap_factor;
gain_rate_eqn.N_total = N_total;
gain_rate_eqn.ASE_factor = ASE_factor;
gain_rate_eqn.E_photon = E_photon;

end