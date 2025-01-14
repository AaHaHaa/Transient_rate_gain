function [N,optimal_dN,...
          extended_Aw,extended_Power_pump] = find_periodic_transient_N(direction,...
                                                                       sim,gain_rate_eqn,...
                                                                       Nt,dt,num_windows,...
                                                                       N0,dN,...
                                                                       At_forward,At_backward,...
                                                                       Power_pump_forward,Power_pump_backward)
%FIND_PERIODIC_TRANSIENT_N It finds the population (under periodic 
%operation of the laser) of N(t)
%
% It finds the populations through MATLAB's optimization function, fminsearch.
% I've initially implemented one that finds the populations by iterating
% between two time windows. However, I find that it converges too slowly.
% In particular, in the close-to-steady-state situations, N(t) varies
% slowly, so iterations between two time widnows also give a slowly varying
% population that leads to slow convergence (to satisfy the periodic
% condition up to a small tolerance value).
% An optimization function is much faster here, which can reduce the
% tolerance 10x smaller after only 10-20 iterations. On the contrary,
% the method with iterations between two windows cannot reduce to half even
% up to 300-500 iterations from our experience/tests.
%
% Input:
%   direction: 'forward' or 'backward';
%              This argument doesn't affect the computation of N(t). It's
%              used to determine which Aw and Power_pump should be output
%              after temporal stretching with dummy zeros for acyclic DFT.
%              In principle, this can be done outside this function, but it
%              saves some computation by outputing the generated parameters
%              in this function anyway.
%   sim: container for simulation parameters
%   gain_rate_eqn: container for rate-eqn parameters
%   Nt: a (1,2) array; the number of points for two windows
%   dt: a (1,2) array; the temporal sampling period (ps)
%   N: a (1,1,1,2) cell array, each containing population at its time window
%   At_forward: a (1,2) cell array, each containing the forward field at its time window
%   At_backward: a (1,2) cell array, each containing the backward field at its time window
%   Power_pump_forward: a (1,1,1,2) cell array, each containing the forward pump power at its time window
%   Power_pump_backward: a (1,1,1,2) cell array, each containing the backward pump power at its time window
%
% Output:
%   N: computed populations
%   extended_Aw: temporally-extended field in the frequency domain
%   extended_Power_pump: temporally-extended pump power (in the time domain)

% Make it with twice the time window to apply the Fourier Transform without aliasing
extended_Aw_forward = cellfun(@(x) ifft(x,gain_rate_eqn.acyclic_conv_stretch(length(x)),1), At_forward, 'UniformOutput',false);
extended_Aw_backward = cellfun(@(x) ifft(x,gain_rate_eqn.acyclic_conv_stretch(length(x)),1), At_backward, 'UniformOutput',false);

extended_Power_pump_forward = cellfun(@(x) gain_rate_eqn.acyclic_zero_padding(x), Power_pump_forward, 'UniformOutput',false);
extended_Power_pump_backward = cellfun(@(x) gain_rate_eqn.acyclic_zero_padding(x), Power_pump_backward, 'UniformOutput',false);

% Output for other uses
switch direction
    case 'forward'
        extended_Aw = extended_Aw_forward;
        extended_Power_pump = extended_Power_pump_forward;
    case 'backward'
        extended_Aw = extended_Aw_backward;
        extended_Power_pump = extended_Power_pump_backward;
end

% =========================================================================
% Find the population N(t) satisfying the periodic boundary condition with
% MATLAB's optimization function, fminsearch
% =========================================================================
% Optimization tolerance
rough_TolFun = gain_rate_eqn.N.periodic.tolf; % for rough optimization
fine_TolFun = rough_TolFun*0.8; % for fine optimization
% Reduce the tolerance if the field is weak
% This step is super important as weak field cannot lead to a consistent
% population result. In particular, during the initial ASE growth, varying
% weak ASE can create different population results, which in turn affects
% its amplification and the final ASE strength. Without a smaller
% tolerance, ASE will fluctuate from iteration to iteration significantly
% (which is tested by me).
weak_field_yes_forward = true;
weak_field_yes_backward = true;
for window_i = 1:num_windows
    if weak_field_yes_forward % if all forward fields are weak
        weak_field_yes_forward = ( weak_field_yes_forward && all(abs(At_forward{window_i}).^2<1e-3) );
    end
    if weak_field_yes_backward % if all backward fields are weak
        weak_field_yes_backward = ( weak_field_yes_backward && all(abs(At_backward{window_i}).^2<1e-3) );
    end
end
if gain_rate_eqn.ignore_ASE % there is no backward field, so don't consider fine searching based on it
    weak_field_yes_backward = false;
end
if weak_field_yes_forward || weak_field_yes_backward % if there is a weak field in either forward or backward direction
    fine_TolFun = fine_TolFun/100;
end

% Because the stopping criteria of Nelder-Mead simplex algorithm used in
% MATLAB's fminsearch rely on the absolute variations of the vectors and
% the function evaluations, rather than the relative difference, it's
% important to normalize/scale them so that the optimization can be
% terminated correctly only when it's close to an optimum, rather than
% simply having values that are too small.
% During optimization for populations, the variations are always small in
% their absolute values, so without correct normalization, optimization
% scheme will almost always terminate at the first evaluation.
if sim.gpu_yes % fminsearch accepts only non-gpuArray inputs
    dN = gather(dN);
end
if ~any(dN) % dN is all zeros
    max_dN = 1; % abort scaling if dN is all zeros

    if length(gain_rate_eqn.energy_levels) > 2 % two-level system doesn't use rough optimization, so it doesn't need "idx"
        idx = N0{1}(1,:)/sum(N0{1}(1,:),2)>0.01; % run rough optimizatin only for those that dominate
    end
else
    max_dN = max(abs(dN(:)));
    dN = dN/max_dN; % normalize it (will be re-scaled back later for optimal_dN)

    if length(gain_rate_eqn.energy_levels) > 2 % two-level system doesn't use rough optimization, so it doesn't need "idx"
        idx = abs(dN)>0.01; % run rough optimizatin only for those that dominate
    end
end

% Small variations to the input guess for the Nelder-Mead simplex algorithm
% in fminsearch.
% If the value isn't zero, vary by 5% first as initial search
% If the value is zero, vary by (an absolute-valued) 1e-6 as initial search
usual_delta = 0.05;
zero_term_delta = 1e-6;

% Use NRMSE0 to set the convergence tolerance as a number smaller than this
% initial value by a preset factor to ensure improvement of finding better
% N(t) satisfying the periodic boundary condition:
%    tolerance=NRMSE0/gain_rate_eqn.N.periodic.tol_ratio
NRMSE0 = calc_N(sim,gain_rate_eqn,...
                Nt,dt,num_windows,...
                N0,...
                extended_Aw_forward,extended_Aw_backward,...
                extended_Power_pump_forward,extended_Power_pump_backward,...
                zeros(length(gain_rate_eqn.energy_levels)-1,1),true(1,length(gain_rate_eqn.energy_levels)-1),max_dN,1);
if ~any(dN) % the steady-state approximation isn't good enough, so NRMSE0 is huge, making the latter tolerance scaling ineffective. We then reduce the NRMSE0 more to aggressively force the optimization to stop at smaller tolerance.
    NRMSE0 = NRMSE0/1e3;
end
if NRMSE0 < 1e-9 % this typically happens for the first z-step when the steady-state approximation is a perfect solution, which happens when the field is too weak
    NRMSE0 = 1;
end

% Rough optimization with only the dominant dN
finished = false;
if length(gain_rate_eqn.energy_levels) > 2 % don't run it for a two-level system; just run fine optimization
    rough_opt_func = @(dN) calc_N(sim,gain_rate_eqn,...
                                  Nt,dt,num_windows,...
                                  N0,...
                                  extended_Aw_forward,extended_Aw_backward,...
                                  extended_Power_pump_forward,extended_Power_pump_backward,...
                                  dN,idx,max_dN,1/NRMSE0);
    option = optimset('TolFun',rough_TolFun);
    %option = optimset('PlotFcns',@optimplotfval,'TolFun',rough_TolFun); % plot the process of optimization
    [optimal_dN,rough_avg_NRMSE] = myfminsearch(rough_opt_func,dN(idx),option,usual_delta,zero_term_delta);
    dN(idx) = optimal_dN;

    usual_delta = 0.005;

    % If rough optimization is pretty good, then there's no need to run the fine optimization.
    if rough_avg_NRMSE < fine_TolFun
        finished = true;

        [~,N] = rough_opt_func(optimal_dN); % calculate N(t) with the optimal N(t=0 in window 1)
        optimal_dN = dN;
    end
end

% Fine optimization with all the dN
if ~finished
    fine_opt_func = @(dN) calc_N(sim,gain_rate_eqn,...
                                 Nt,dt,num_windows,...
                                 N0,...
                                 extended_Aw_forward,extended_Aw_backward,...
                                 extended_Power_pump_forward,extended_Power_pump_backward,...
                                 dN,true(length(gain_rate_eqn.energy_levels)-1,1),max_dN,1/NRMSE0);
    option = optimset('TolFun',fine_TolFun);
    %option = optimset('PlotFcns',@optimplotfval,'TolFun',fine_TolFun); % plot the process of optimization
    [optimal_dN,fine_avg_NRMSE] = myfminsearch(fine_opt_func,dN,option,usual_delta,zero_term_delta); %#ok

    [~,N] = fine_opt_func(optimal_dN); % calculate N(t) with the optimal N(t=0 in window 1)
end

optimal_dN = optimal_dN*max_dN; % scale it back to population's unit

end

%% Optimization function to find N(t)
function [avg_NRMSE,N] = calc_N(sim,gain_rate_eqn,...
                                Nt,dt,num_windows,...
                                N,...
                                extended_Aw_forward,extended_Aw_backward,...
                                extended_Power_pump_forward,extended_Power_pump_backward,...
                                dN,...
                                idx,...
                                dN_scaling,...
                                NRMSE_scaling) % initial guess for population at t=0 in window 1)
%CALC_N It calculates the population N(t) with the initial N(t=0 in window 1)
% and outputs the error between N(t=0 in window 1) and N(t=t_end in window 2),
% as well as the computed population.
%
% Note the gap between time windows numerically.
% Each time window contains Nt points, which represents the time 
% t = 0, 1, 2,..., (Nt-1) time-unit
% However, the next time window starts with t=Nt, so we need to extrapolate
% to obtain the value at t=Nt for the continuity.

N0 = N{1}(1,:);
N0(:,idx) = N0(:,idx) + dN.'*dN_scaling*gain_rate_eqn.N_total;
if sum(N0) > gain_rate_eqn.N_total
    N0 = N0/sum(N0)*gain_rate_eqn.N_total;
end

% ---------------------------------------------------------------------
% Without coherent pulses (time window 1)
% ---------------------------------------------------------------------
N{1} = N_scaling(N{1}(1,:),N0,N{1},gain_rate_eqn.N_total);
N{1} = solve_gain_rate_eqn(sim,gain_rate_eqn,...
                           gain_rate_eqn.acyclic_zero_padding(N{1}),...
                           extended_Aw_forward{1},extended_Aw_backward{1},...
                           extended_Power_pump_forward{1},extended_Power_pump_backward{1},...
                           gain_rate_eqn.extended_cross_sections{1},gain_rate_eqn.extended_E_photon{1},...
                           Nt,dt(1));
% ---------------------------------------------------------------------
% With and without coherent pulses (other time windows)
% ---------------------------------------------------------------------
for window_i = 2:num_windows
    NpreviousEnd = interp1((0:Nt-1)',N{window_i-1},Nt,'linear','extrap');
    N{window_i} = N_scaling(N{window_i}(1,:),NpreviousEnd,N{window_i},gain_rate_eqn.N_total);
    N{window_i} = solve_gain_rate_eqn(sim,gain_rate_eqn,...
                                      gain_rate_eqn.acyclic_zero_padding(N{window_i}),...
                                      extended_Aw_forward{window_i},extended_Aw_backward{window_i},...
                                      extended_Power_pump_forward{window_i},extended_Power_pump_backward{window_i},...
                                      gain_rate_eqn.extended_cross_sections{mod(window_i+1,2)+1},gain_rate_eqn.extended_E_photon{mod(window_i+1,2)+1},...
                                      Nt,dt(window_i));
end

weight = N0/sum(N0,2); weight(weight<1e-3)=0;
NfinalEnd = interp1((0:Nt-1)',N{end},Nt,'linear','extrap');
NRMSE = abs(N0 - NfinalEnd)./((N0+NfinalEnd)/2).*weight;
NRMSE(isnan(NRMSE)) = 0; % exclude levels with zero population
avg_NRMSE = sum(NRMSE,2)*NRMSE_scaling;

if sim.gpu_yes
    avg_NRMSE = gather(avg_NRMSE);
end

end

%% Helper function
function y = N_scaling(x0,x1,y0,x_total)
%N_SCALING It re-scales populations with the population at t_end of the
%other time window with an exponential function.
%
% I have initially implemented a simple replacement by replacing the first
% element only through these two lines after N{2} and N{1} computations:
%   N{1}(1,:) = N{2}(end,:)
%   and
%   N{2}(1,:) = N{1}(1,:)
% However, this simple replacement creates a discrete jump between the
% first and the second elements in time, potentially resulting in difficult
% convergence of MPA. Therefore, I implement this re-scaling function to
% scale the new N1(all t) corresponding to the N2(t_end) and N1(t_first),
% or vice versa.
%
%   x0: N1 at its t_first
%   x1: N2 at its t_end
%   y0: the entire population [ x0=y0(t_first) ]
%   x_total: total population
%   y: scaled population
%
% The monotonically-scaling function follows 
%   if y0(t) > x0, it's scaled exponentially from x1 to 1 based on the relation between x0 and x1
%   if y0(t) < x0, it's scaled exponentially from 0 to x1
%   if y0(t) = x0, it becomes x1

% Make them population ratio with 0-1 for scaling
x0 = x0/x_total;
x1 = x1/x_total;
y0 = y0/x_total;

c = 6; % a coefficient controlling the slope of scaling; just a random number I found to be useful
y = y0;
for Nidx = 1:size(x0,2)
    idx = (y0(:,Nidx) >= x0(:,Nidx));
    
    y( idx,Nidx) = (1-x1(:,Nidx)).*(1 - exp(-c*(y0( idx,Nidx)-x0(:,Nidx)))) + x1(:,Nidx);
    y(~idx,Nidx) =   -x1(:,Nidx) .*(1 - exp( c*(y0(~idx,Nidx)-x0(:,Nidx)))) + x1(:,Nidx);
end

% Ensure that (sum(y,2) < total population) after scaling
sum_y = sum(y,2);
idx = sum_y>1;
y(idx,:) = y(idx,:)./sum_y(idx);

% Make it into population, rather than ratio
y = y*x_total;

end