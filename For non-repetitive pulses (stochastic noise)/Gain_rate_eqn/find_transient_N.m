function [N,...
          extended_Aw_forward,extended_Power_pump_forward] = find_transient_N(sim,gain_rate_eqn,...
                                                                              Nt,dt,num_windows,...
                                                                              N,...
                                                                              At_forward,...
                                                                              Power_pump_forward)
%FIND_TRANSIENT_N It finds the population of N(t)
%
% Input:
%   sim: container for simulation parameters
%   gain_rate_eqn: container for rate-eqn parameters
%   Nt: a (1,1,1,num_windows) array; the number of points for two windows
%   dt: a (1,1,1,num_windows) array; the temporal sampling period (ps)
%   num_windows: the number of time windows
%   N: a (1,1,1,num_windows) cell array, each containing population at its time window
%   At_forward: a (1,1,1,num_windows) cell array, each containing the forward field at its time window
%   Power_pump_forward: a (1,1,1,num_windows) cell array, each containing the forward pump power at its time window
%
% Output:
%   N: computed populations
%   extended_Aw_forward: temporally-extended field in the frequency domain
%   extended_Power_pump_forward: temporally-extended pump power (in the time domain)

% Make it with twice the time window to apply the Fourier transform without aliasing
extended_Aw_forward = cellfun(@(x) ifft(x,gain_rate_eqn.acyclic_conv_stretch(length(x)),1), At_forward, 'UniformOutput',false);

extended_Power_pump_forward = cellfun(@(x) gain_rate_eqn.acyclic_zero_padding(x), Power_pump_forward, 'UniformOutput',false);

% -------------------------------------------------------------------------
% Compute photon rate, which is fixed during searching for periodic N(t)
% -------------------------------------------------------------------------
h = 6.62607004e-34;
c = 299792458;

R = cell(1,1,1,num_windows);
for window_i = 1:num_windows
    % Compute F^(-1)[A(Omega) * sqrt(sigma/hbar/omega)]
    At_forward_i = fft(extended_Aw_forward{window_i}.*sqrt(gain_rate_eqn.extended_cross_sections{mod(window_i+1,2)+1}./gain_rate_eqn.extended_E_photon{mod(window_i+1,2)+1}),[],1); % sqrt(m^2/s); size: (Nt,num_cross_sections)
    
    A2_i = abs(At_forward_i).^2; % A2(t); the number of scattered photons

    % Rate/photon_energy:
    % size: (Nt,num_cross_sections), unit: W
    R_pump  = gain_rate_eqn.overlap_factor.pump.*gain_rate_eqn.cross_sections_pump.*extended_Power_pump_forward{window_i} /h*(gain_rate_eqn.pump_wavelength*1e-9)/c; % 1/s
    
    R_ASE_signal = gain_rate_eqn.overlap_factor.signal.*A2_i; % 1/s; size: (Nt,num_cross_sections)
    
    % Total photon rate
    R{window_i} = R_pump + R_ASE_signal;
end

% ---------------------------------------------------------------------
% Without coherent pulses (time window 1)
% ---------------------------------------------------------------------
N{1} = solve_gain_rate_eqn(sim,gain_rate_eqn,...
                           gain_rate_eqn.acyclic_zero_padding(N{1}),...
                           R{1},...
                           Nt,dt(1));

% ---------------------------------------------------------------------
% With and without coherent pulses (other time windows)
% ---------------------------------------------------------------------
for window_i = 2:num_windows
    NpreviousEnd = interp1((0:Nt-1)',N{window_i-1},Nt,'linear','extrap');
    N{window_i} = N_scaling(N{window_i}(1,:),NpreviousEnd,N{window_i},gain_rate_eqn.N_total);
    N{window_i} = solve_gain_rate_eqn(sim,gain_rate_eqn,...
                                      gain_rate_eqn.acyclic_zero_padding(N{window_i}),...
                                      R{window_i},...
                                      Nt,dt(window_i));
end

end

%% Helper function
function y = N_scaling(x0,x1,y0,x_total)
%N_SCALING It re-scales populations with the population at t_end of the
%other time window with an exponential function.
%
% I have initially implemented a simple replacement by replacing the first
% element only through these two lines after N{2} and N{1} computations for
% the initial two-time-window demonstration:
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