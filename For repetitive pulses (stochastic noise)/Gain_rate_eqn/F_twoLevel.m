function dN1dt = F_twoLevel(N1,N_total,A,Gamma,kijkl,R) %#ok; suppress warning for unused variables
%F_TWOLEVEL container for the coupled equations of population evolutions
%for a two-level system
%
% Input arguments:
%   N0: ground-state population; (Nt,1)
%   N_total: total population; a scalar
%   A: 1/tau in this two-level-system function, tau: lifetime of the excited state.
%      This is the summation of spontaneous transition rate/probabilities and multiphonon decay rate; a scalar
%   Gamma: multiphonon decay rates (unused here)
%   kijkl: nonlinear coupling coefficients for upconversion or cross relaxation effects (unused here)
%   R: photon rate; (Nt,num_cross_sections), num_cross_sections = 2
%
% Output argument:
%   dN1dt: variation of populations; (Nt,1)

% Populations of each energy level
%N0 = N_total - N1; % ground state

% Stimulated photon rates (absorption and emission)
wa01 = R(:,1);
we10 = R(:,2);

% Variation of populations
% dN0 + dN1 = dN_total = 0
%dN0dt = A*N1 ...
%        + we10.*N1 ...
%        - wa01.*N0;
%dN1dt = -A*N1 ...
%        - we10.*N1 ...
%        + wa01.*N0;
dN1dt = -A*N1 ...
        - (we10+wa01).*N1 ...
        + wa01.*N_total;

end