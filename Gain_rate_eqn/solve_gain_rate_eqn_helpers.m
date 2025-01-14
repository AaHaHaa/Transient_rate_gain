function func = solve_gain_rate_eqn_helpers
%SOLVE_GAIN_RATE_EQN_HELPERS containers for the rate-equation solvers

func.solve_N = @solve_N; % solve the population
func.solve_pump = @solve_pump; % solve the pump
func.solve_Power = @solve_Power; % solve the pump and ASE powers, and the signal gain

func.solve_N_SS = @solve_N_SS; % solve the steady-state population as an initial guess for the algorithm

end

%% solve_N
function N = solve_N(sim,gain_rate_eqn,...
                     Nt,dt,...
                     N,N_total,...
                     overlap_factor,cross_sections_pump,...
                     A2,...
                     Power_pump_forward,Power_pump_backward)
%SOLVE_N It solves the populations among various energy levels.
%
% computational dimension here: (Nt,num_cross_sections/num_levels)

h = 6.62607004e-34;
c = 299792458;

% Rate/photon_energy:
% size: (Nt,num_cross_sections), unit: W
num_cross_sections = length(cross_sections_pump);
if sim.gpu_yes
    R_forward  = zeros([Nt,num_cross_sections],'gpuArray');
    R_backward = zeros([Nt,num_cross_sections],'gpuArray');
else
    R_forward  = zeros([Nt,num_cross_sections]);
    R_backward = zeros([Nt,num_cross_sections]);
end
if ismember(gain_rate_eqn.pump_direction,{'co','bi'})
    R_forward  = overlap_factor.pump.*cross_sections_pump.*Power_pump_forward /h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
end
if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
    R_backward = overlap_factor.pump.*cross_sections_pump.*Power_pump_backward/h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
end
R.pump = R_forward + R_backward;

R.ASE_signal = overlap_factor.signal.*A2; % 1/s; size: (Nt,num_cross_sections)

% Total photon rate
R = R.pump + R.ASE_signal;

% Compute the populations with the massively parallel algorithm (MPA).
% The derivative dN/dt is found at all time points simultaneously with the 
% previously known N(t), which is later used to compute the next N(t) for 
% all time points at once. With the next N(t), a newer dN/dt at all time 
% points is found for the next N(t). This process continues until N(t) at 
% all time points converge to steady-state values.
% However, if N(t) cannot converge within preset maximum iterations
% (= gain_rate_eqn.N_iter.max_iter), N(t) is computed time point by time
% point, rather than MPA:
%   N(t+dt) = N(t) + ( dN/dt(t) )*dt
% This method doesn't require any convergence mechanism as MPA. The reason
% MPA is preferred because it's much faster. Typical iterations for MPA
% convergence is as small as 3. I found that the time-point-by-time-point
% algorithm is triggered only when
% 1. incoherent ASE is so strong such that N(t) and dN/dt(t) are highly
%    fluctuating, making MPA difficult to converge.
% 2. Time window is too huge, compared to the characteristic times of
%    energy levels. Since MPA can only converge when variations are small
%    within the time window considered, it fails if the time window is too 
%    large. I also found this to happen more often for multi-level systems
%    even with a smaller time window. The complciated governing equation
%    doesn't ensure MPA convergence.
F = @(N,R) gain_rate_eqn.N.eqn.tr(N,N_total,gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.N.eqn.kijkl,R);
last_N = N;
avg_NRMSE = zeros(1,gain_rate_eqn.N.MPA.max_iter);
for n_iter = 1:gain_rate_eqn.N.MPA.max_iter % maximum iterations
    % ---------------------------------------------------------------------
    % Compute N(t) through massively parallel algorithm (MPA)
    % ---------------------------------------------------------------------
    dNdt = F(N,R);
    dNdt(ceil(Nt/2)+1:end,:) = 0; % remove unused information in the extended time window
    
    % Integrate to obtain N(t) from dNdt
    dNdt(1,:) = dNdt(1,:)/2;
    N(2:end,:) = N(1,:) + (cumsum(dNdt(1:end-1,:),1) + dNdt(2:end,:)/2)*(dt*1e-12);

    % Set boundaries
    N(ceil(Nt/2)+1:end,:) = 0; % remove unused information in the extended time window
    N(N<0) = 0; % lower bound: population cannot be negative
    idx = sum(N,2)>N_total;
    N(idx,:) = N(idx,:)./sum(N(idx,:),2)*N_total; % upper bound: in case it blows up during MPA comptutations, which can easily happen for multi-level computations

    % Error computation with normalized root-mean-square error
    NRMSE = sqrt(mean(abs(N-last_N).^2,1))./N_total;
    NRMSE(isnan(NRMSE)) = 0; % exclude where there is no doping, N_total=0
    avg_NRMSE(n_iter) = sum(NRMSE);
    
    % If it has converged to within tol, then quit
    if n_iter > 1
        avg_NRMSE_convergence_rate = abs(avg_NRMSE(n_iter)/avg_NRMSE(n_iter-1)-1);
        if ((avg_NRMSE_convergence_rate < gain_rate_eqn.N.MPA.tol && avg_NRMSE(n_iter) < 1e-3) || ... % population variation is not only small but also converging
            avg_NRMSE(n_iter) == 0) && ... % population convergence rate becomes zero
           n_iter >= gain_rate_eqn.N.MPA.min_iter % minimum iterations
            break
        end
    end


    % ---------------------------------------------------------------------
    % If MPA fails to converge, then
    % compute N(t) through sequential integration (with RK4): N(t+dt) = N(t) + ( dN/dt(t) )*dt
    % ---------------------------------------------------------------------
    if n_iter > gain_rate_eqn.N.MPA.max_iter/gain_rate_eqn.N.MPA.start_sequential_factor
        if avg_NRMSE(n_iter) > 0.1 || ... % avg_NRMSE is just too huge so that we don't expect MPA to converge; most likely, it is diverging
           n_iter == gain_rate_eqn.N.MPA.max_iter % maximum iterations
            % N is solved one by one, so it's slower (by a lot) with GPU due to data creation and moving overheads
            if sim.gpu_yes
                N = gather(N);
                R = gather(R);
            end
            for ti = 2:ceil(Nt/2)
                % Solve with RK4 (explicit fourth-order Runge-Kutta)
                a1 = F(N(ti-1,:)                  ,R(ti-1,:));
                a2 = F(N(ti-1,:)+a1*((dt*1e-12)/2),R(ti-1,:));
                a3 = F(N(ti-1,:)+a2*((dt*1e-12)/2),R(ti-1,:));
                a4 = F(N(ti-1,:)+a3*((dt*1e-12))  ,R(ti-1,:));
                N(ti,:) = N(ti-1,:) + (a1+2*a2+2*a3+a4)*((dt*1e-12)/6);
                % Solve with the Trapezoidal method (implicit quadrature)
                %dNdt = F(N(ti-1,:),R(ti-1,:));
                %N(ti,:) = N(ti-1,:) + dNdt*(dt*1e-12);
    
                % Set boundaries
                N(ti,N(ti,:)<0) = 0; % lower bound: population cannot be negative
                if sum(N(ti,:),2) > N_total
                    N(ti,:) = N(ti,:)/sum(N(ti,:),2)*N_total; % upper bound: total population is N_total
                end
            end
            break; % if force_sequential=true, this block of code will run when n_iter=1 and finishes computing N
        end
    end
    % ---------------------------------------------------------------------

    last_N = N; % transfer data to last_N for the next MPA convergence check
end

end

%% %% Contained function: solve_Power
function dPdz_pump_forward = solve_pump(gain_rate_eqn,...
                                        N,...
                                        P0)

num_windows = length(P0);
dPdz_pump_forward = cell(1,1,1,num_windows);
for window_i = 1:num_windows
    dPdz_pump_forward{window_i} = solve_Power('pump',...
                                              gain_rate_eqn,...
                                              [],...
                                              gain_rate_eqn.cross_sections_pump,...
                                              N{window_i},...
                                              P0{window_i},[],[]);
end

end

%% Helper function: solve_Power
function varargout = solve_Power(field_type,...
                                 gain_rate_eqn,...
                                 dz,...
                                 cross_sections,...
                                 N,...
                                 P0,SE_power,Aw)
%SOLVE_POWER solves Power(z+dz) for pump and ASE+signal.
%
%   dz: um
%

Nt = size(N,1);
N = cat(2,gain_rate_eqn.N_total-sum(N,2),N); % include the gound-state population
N(gain_rate_eqn.acyclic_conv_shrink(Nt)+1:end,:) = 0; % remove population in the extended time window

A_core = pi*(gain_rate_eqn.core_diameter/2)^2;
switch field_type
    case {'signal','sponASE'}
        overlap_factor = gain_rate_eqn.overlap_factor.signal;
    case 'pump'
        overlap_factor = gain_rate_eqn.overlap_factor.pump;
end
overlap_factor = overlap_factor*A_core; % For the single mode, the integral w.r.t. x and y can be done first, which becomes overlap_factor here.

switch field_type
    case 'signal'
        gAt = overlap_factor*sum(gain_rate_eqn.plusminus.*N(:,gain_rate_eqn.N_idx).*fft(cross_sections.*Aw,[],1),2)/2*1e6;
        
        varargout = {gAt};
    case 'sponASE'
        spon_idx = (squeeze(gain_rate_eqn.plusminus)>0);
        dAt_SE = sqrt(overlap_factor)*sum(gain_rate_eqn.plusminus(spon_idx).*sqrt(N(:,gain_rate_eqn.N_idx(spon_idx)))...
                                          .*fft(sqrt(cross_sections(:,spon_idx).*SE_power*1e6*dz*gain_rate_eqn.ASE_factor).*randn(Nt,sum(spon_idx)).*exp(1i*2*pi*rand(Nt,sum(spon_idx)))),2);
        
        varargout = {dAt_SE};
    case 'pump'
        dPdz = overlap_factor*sum(gain_rate_eqn.plusminus.*cross_sections.*N(:,gain_rate_eqn.N_idx),2).*P0*1e6; % unit: W
        
        varargout = {dPdz};
end

end

%% ========================================================================
% Below are steady-state computation functions for finding the initial guess
% =========================================================================

%% Contained function: solve_N
function N = solve_N_SS(gain_rate_eqn,gpu_yes,...
                        N_total,...
                        E_photon,...
                        overlap_factor,cross_sections_pump,cross_sections,...
                        Power_pump,...
                        A2)
%SOLVE_N_SS It solves the steady-state populations among various energy
%levels.
%
% This function is a container of the population solver.
% The solver is split into two helper functions, one for only-two-level
% systems and the other for the general multi-level systems. Because
% two-level systems don't have complicated nonlinear terms, such as
% cross-relaxation effects in Er, Tm, etc., their population can be easily
% solved with
%   N2 = absorption_term / (absorption_term + emission_term + 1/tau), tau: upper-state lifetime
%
% However, for multi-level systems, a numerical solver is required. I
% implemented the trust-region method with a dogleg sub-problem solver.
%
% Due to these two processes, it's helpful to split the population solver 
% into two functions based on whether it's a two-level system or not.

if length(gain_rate_eqn.energy_levels) == 2 % two-level system
    solve_N_func = str2func('solve_N_SS_2levels');
else % multi-level system
    solve_N_func = str2func('solve_N_SS_levels');
end

N = solve_N_func(gain_rate_eqn,gpu_yes,...
                 N_total,...
                 E_photon,...
                 overlap_factor,cross_sections_pump,cross_sections,...
                 Power_pump,...
                 A2);

end

function N1 = solve_N_SS_2levels(gain_rate_eqn,gpu_yes,...
                                 N_total,...
                                 E_photon,...
                                 overlap_factor,cross_sections_pump,cross_sections,...
                                 Power_pump,...
                                 A2)
%SOLVE_N_SS_2LEVELS solves the ion density in the upper state under the 
%steady state of the rate equation

h = 6.62607004e-34;
c = 299792458;

% Rate/photon_energy:
R_over_photon.pump  = overlap_factor.pump.*cross_sections_pump.*Power_pump /h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s

% It's real! Use "real" to save the memory.
R_over_photon.signal = real(overlap_factor.signal*sum(A2.*(cross_sections./E_photon),1)); % 1/s; size: (1,[GSA01,emi10])

% ion density in the upper state
% N2 = N_total*sum(all the absorption terms)/
%      ( sum(all the absorption and emission terms) + upper-state-lifetime term )
total_absorption = R_over_photon.pump(:,1) + R_over_photon.signal(:,1);
total_emission   = R_over_photon.pump(:,2) + R_over_photon.signal(:,2);
N1 = N_total.*total_absorption./... % size: (Nx,Nx,1,1,1,M)
     (total_absorption + total_emission + ... % absorption and emission terms
      1/gain_rate_eqn.N.eqn.tau);                   % upper-state-lifetime term

end

function N = solve_N_SS_levels(gain_rate_eqn,gpu_yes,...
                               N_total,...
                               E_photon,...
                               overlap_factor,cross_sections_pump,cross_sections,...
                               Power_pump,...
                               A2)
%SOLVE_N_SS_LEVELS It solves the steady-state populations among various 
%energy levels.

h = 6.62607004e-34;
c = 299792458;

% Rate/photon_energy:
R_over_photon.pump  = overlap_factor.pump.*cross_sections_pump.*Power_pump /h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s

% It's real! Use "real" to save the memory.
R_over_photon.signal = real(overlap_factor.signal*sum(A2.*(cross_sections./E_photon),1)); % 1/s; size: (1,num_cross_sections)

% Compute the steady-state populations
R_over_photon = shiftdim(R_over_photon.pump + R_over_photon.signal,-8); % shift the dimension to prepare for the solver
N = [gain_rate_eqn.N_total;zeros(length(gain_rate_eqn.energy_levels)-2,1)];
F = @(N) gain_rate_eqn.N.eqn.ss.F(N,N_total,gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.N.eqn.kijkl,R_over_photon);
J = @(N) gain_rate_eqn.N.eqn.ss.J(N,N_total,gain_rate_eqn.N.eqn.Arad,gain_rate_eqn.N.eqn.Gammai,gain_rate_eqn.N.eqn.kijkl,R_over_photon);
N = myTrustRegion(F,J,N,10,1e-6,size(R_over_photon,8),N_total,gpu_yes);

end