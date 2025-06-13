function [A1t_forward,...
          Power_pump_forward,...
          N,dN] = forward_stepping_RK4IP_rategain(A0t_forward,A0t_backward,...
                                                  Power_pump_forward,Power_pump_backward,...
                                                  N,dN,...
                                                  Nt,dt,num_windows,...
                                                  sim,gain_rate_eqn,...
                                                  SK_info,SRa_info,SRb_info,...
                                                  haw,hbw,...
                                                  n2_prefactor,expD,...
                                                  At_noise)
%FORWARD_STEPPING_RK4IP_RATEGAIN Take one step with RK4IP with a gain model  
%solved from rate equations.
%
% Input:
%    A0t_forward - initial forward-propagating field in the time domain
%                  a (1,1,1,num_windows) cell array, each containing a (Nt,1) array; sqrt(W)
%    A0t_backward - initial backward-propagating field in the time domain
%                  a (1,1,1,num_windows) cell array, each containing a (Nt,1) array; sqrt(W)
%
%    Power_pump_forward - the power of the co-propagating pump
%                         a (1,1,1,num_windows) cell array, each containing a (Nt,1) array; in W
%    Power_pump_backward - the power of the counter-propagating pump
%                          a (1,1,1,num_windows) cell array, each containing a (Nt,1) array; in W
%
%    N - populations
%        a (1,1,1,num_windows) cell array, each containing a (Nt,num_levels-1) array; in m^(-3)
%    dN - (1,num_levels-1); the population variations over one dz
%
%    Nt - the number of sampling points
%    dt - time grid point spacing; ps;  a (1,num_windows) array
%    num_windows - the number of windows
%
%    sim - container for simulation parameters
%    gain_rate_eqn - container for rate-eqn parameters
%
%    SRa_info.SRa - SRa tensor; m^-2
%    SRa_info.nonzero_midx1234s - required SRa indices in total
%    SRa_info.nonzero_midx34s - required (SRa) indices for partial Raman term (only for CPU computation)
%    SRb_info.SRb - SRb tensor; m^-2
%    SRb_info.nonzero_midx1234s - required SRb indices in total
%    SRb_info.nonzero_midx34s - required (SRb) indices for partial Raman term (only for CPU computation)
%    SK_info.SK - SK tensor; m^2 (unempty if considering polarizaton modes)
%    SK_info.nonzero_midx1234s - required SK indices in total (unempty if considering polarizaton modes)
%
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    n2_prefactor - the nonlinear prefactor = 1i*n2*omega/c; m/W; (Nt,1)
%    expD - dispersion term exp(i*D*dz/2); (Nt,1)
%
%    At_noise - the shot noise in the windows
%               a (1,num_windows) cell array, each containing a (Nt,1) array; sqrt(W)
%
% Output:
%    A1t_forward - the forward-propagating field in the time domain after one step size
%                  a (1,1,1,num_windows) cell array, each containing a (Nt,1) array; sqrt(W)
%    Power_pump_forward - the power of the co-propagating pump
%                         a (1,1,1,num_windows) cell array, each containing a (Nt,1) array; in W
%    N - populations
%        a (1,1,1,num_windows) cell array, each containing a (Nt,num_levels-1) array; in m^(-3)
%    dN - (1,num_levels-1); the population variations over one dz
                            
A1t_forward = cell(1,1,1,num_windows);

% =========================================================================
% Find the population N(t) under the periodic boundary condition
% =========================================================================
[N,dN,...
 extended_A0w_forward,extended_Power_pump_forward] = find_periodic_transient_N('forward',...
                                                                               sim,gain_rate_eqn,...
                                                                               Nt,dt,num_windows,...
                                                                               N,dN,...
                                                                               A0t_forward,A0t_backward,...
                                                                               Power_pump_forward,Power_pump_backward);

% =========================================================================
% Finished finding the population N and start to compute the field (signal 
% fields and ASE) and the pump power
% =========================================================================

extended_N = cellfun(@(x) gain_rate_eqn.acyclic_zero_padding(x), N,'UniformOutput',false);

func = solve_gain_rate_eqn_helpers();
dPdz_pump_forward = func.solve_pump(gain_rate_eqn,...
                                    extended_N,...
                                    extended_Power_pump_forward);

% -------------------------------------------------------------------------
% Without coherent pulses
% -------------------------------------------------------------------------

for window_i = 1:2:num_windows
    % Pump power is found with a simple Newton's method
    Power_pump_forward{window_i} = Power_pump_forward{window_i} + dPdz_pump_forward{window_i}(1:Nt)*sim.dz;
    % Pump power cannot be negative
    Power_pump_forward{window_i}(Power_pump_forward{window_i}<0) = 0;

    dAt_SE1 = func.solve_Power('sponASE',...
                               gain_rate_eqn,...
                               sim.dz,...
                               gain_rate_eqn.extended_cross_sections{1},...
                               extended_N{window_i},...
                               [],abs(ifft(At_noise{window_i},[],1)).^2,[]);
    % Field is found with RK4
    A1t_forward{window_i} = RK4_ASE(extended_A0w_forward{window_i},...
                                    gain_rate_eqn,...
                                    sim.dz,...
                                    extended_N{window_i},...
                                    dAt_SE1);
end

% -------------------------------------------------------------------------
% With coherent pulses
% -------------------------------------------------------------------------

for window_i = 2:2:num_windows
    % Pump power is found with a simple Newton's method
    Power_pump_forward{window_i} = Power_pump_forward{window_i} + dPdz_pump_forward{window_i}(1:Nt)*sim.dz;
    % Pump power cannot be negative
    Power_pump_forward{window_i}(Power_pump_forward{window_i}<0) = 0;

    dAt_SE2 = func.solve_Power('sponASE',...
                               gain_rate_eqn,...
                               sim.dz,...
                               gain_rate_eqn.extended_cross_sections{2},...
                               extended_N{window_i},...
                               [],abs(ifft(At_noise{window_i},[],1)).^2,[]);

    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    % Field is found with RK4IP
    % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    % Set up matrices for the following Kerr, Ra, and Rb computations
    if sim.gpu_yes
        Kerr = complex(zeros(gain_rate_eqn.acyclic_conv_stretch(Nt),1,1, 'gpuArray'));
        Ra = complex(zeros(gain_rate_eqn.acyclic_conv_stretch(Nt),1,1, 'gpuArray'));
        Rb = complex(zeros(gain_rate_eqn.acyclic_conv_stretch(Nt),1,1, 'gpuArray'));
    else
        Kerr = complex(zeros(gain_rate_eqn.acyclic_conv_stretch(Nt),1));
        Ra = complex(zeros(gain_rate_eqn.acyclic_conv_stretch(Nt),1,1));
        Rb = complex(zeros(gain_rate_eqn.acyclic_conv_stretch(Nt),1,1));
    end

    % Propagate through the first dispersion + gain section
    A_IP = expD.*extended_A0w_forward{window_i};

    % Propagate through the nonlinearity with RK4
    a1 = expD.*N_op(       extended_A0w_forward{window_i},...
                    sim,gain_rate_eqn,...
                    SK_info, SRa_info, SRb_info,...
                    Kerr, Ra, Rb,...
                    haw, hbw,...
                    At_noise{window_i},...
                    n2_prefactor,...
                    gain_rate_eqn.acyclic_conv_stretch(Nt),...
                    extended_N{window_i});
    a2 =       N_op(       A_IP+a1*(sim.dz/2),...
                    sim,gain_rate_eqn,...
                    SK_info, SRa_info, SRb_info,...
                    Kerr, Ra, Rb,...
                    haw, hbw,...
                    At_noise{window_i},...
                    n2_prefactor,...
                    gain_rate_eqn.acyclic_conv_stretch(Nt),...
                    extended_N{window_i});
    a3 =       N_op(       A_IP+a2*(sim.dz/2),...
                    sim,gain_rate_eqn,...
                    SK_info, SRa_info, SRb_info,...
                    Kerr, Ra, Rb,...
                    haw, hbw,...
                    At_noise{window_i},...
                    n2_prefactor,...
                    gain_rate_eqn.acyclic_conv_stretch(Nt),...
                    extended_N{window_i});
    a4 =       N_op(expD.*(A_IP+a3*(sim.dz)),...
                    sim,gain_rate_eqn,...
                    SK_info, SRa_info, SRb_info,...
                    Kerr, Ra, Rb,...
                    haw, hbw,...
                    At_noise{window_i},...
                    n2_prefactor,...
                    gain_rate_eqn.acyclic_conv_stretch(Nt),...
                    extended_N{window_i});

    A1 = expD.*(A_IP + (a1+2*a2+2*a3)*(sim.dz/6)) + a4*(sim.dz/6);

    A1t_forward{window_i} = fft(A1,[],1) + dAt_SE2; % add spontaneous emission
end

% Downsample it back to the original time window
for window_i = 1:num_windows
    A1t_forward{window_i} = A1t_forward{window_i}(1:Nt);
end

end

%% Nonlinear operator
function dAwdz = N_op(Aw,...
                      sim,gain_rate_eqn,...
                      SK_info, SRa_info, SRb_info,...
                      Kerr, Ra, Rb,...
                      haw, hbw,...
                      At_noise,...
                      n2_prefactor,...
                      Nt,...
                      N)
%N_op Calculate dAdz

At = fft(Aw,[],1);
At_wNoise = At + At_noise;

% Calculate large num_modes^4 Kerr, Ra, and Rb terms.
% If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
if sim.gpu_yes
    % If using the GPU, do the computation with fast CUDA code
    if sim.scalar % scalar fields
        [Kerr,...
         Ra] = feval(sim.cuda_SRSK,...
                     Kerr, Ra,...
                     complex(At_wNoise),...
                     SK_info.SK, SRa_info.SRa,...
                     SRa_info.nonzero_midx1234s,...
                     SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                     sim.include_Raman,...
                     int32(Nt), 1,...
                     1,...
                     sim.cuda_num_operations_SRSK);
    else % polarized fields
        [Kerr,...
         Ra, Rb] = feval(sim.cuda_SRSK,...
                         Kerr, Ra, Rb,...
                         complex(At_wNoise),...
                         SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                         SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                         SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                         sim.include_Raman, ~isempty(hbw),...
                         int32(Nt), 1,...
                         1,...
                         sim.cuda_num_operations_SRSK);
    end
    Kerr = sum(Kerr,3);
else
    % If using the CPU, first precompute Ra_mn and Rb_mn.
    if sim.include_Raman
        midx34s_sub2ind = @(x)...
            cellfun(@(xx)...
                feval(@(sub) sub2ind(1*ones(1,2),sub{:}), num2cell(xx)),... % this extra "feval" is to get "xx", which is of the size 2x1, into the input arguments of "sub2ind", so transforming "xx" into a 2x1 cell, each containing an integer, and using {:} expansion is necessary
            mat2cell(x,2,ones(1,size(x,2)))); % transform (2,num_nonzero34) midx34s into linear indices of a num_modes-by-num_modes matrix
            % What "midx34s_sub2ind" does (e.g.):
            %
            %   x = [1 3;
            %        5 4]
            %
            %   After "mat2cell": {[1;  {[3;  (2x1 cells, each having 2x1 array)
            %                       5]}   4]}
            %
            %   First,
            %
            %   xx = {[1;  , then after "num2cell": {{1}; (1 cell with 2x1 cell)
            %          5]}                           {5}}
            %
            %   The purpose of separating 1 and 5 into cells is to use
            %   index expansion, {:}, to put them into the input
            %   arguments of "sub2ind" function.
            %
            %   For 6 modes and thus for 6x6 matrix, sub2ind([6 6],1,5) = 25
            %
            %   Do the same for xx = {[3;  and get sub2ind([6 6],3,4) = 21
            %                          4]}
            %   Finally, midx34s_sub2ind = [25 21] (1x2 array)

        SRa_nonzero_midx34s = midx34s_sub2ind(SRa_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
        Ra_mn = At_wNoise(:, SRa_info.nonzero_midx34s(1,:)).*conj(At_wNoise(:, SRa_info.nonzero_midx34s(2,:))); % (N,num_nonzero34)
        if ~isempty(hbw)
            SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
            Rb_mn = At_wNoise(:, SRb_info.nonzero_midx34s(1,:)).*conj(At_wNoise(:, SRb_info.nonzero_midx34s(2,:))); % (N,num_nonzero34)
        end
    end
    
    % Then calculate Kerr, Ra, and Rb.
    for midx1 = 1
        % Kerr
        nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
        midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
        midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
        midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
        Kerr(:,midx1) = sum(permute(SK_info.SK(nz_midx1),[2 1]).*At_wNoise(:, midx2).*At_wNoise(:, midx3).*conj(At_wNoise(:, midx4)),2);
        if sim.include_Raman
            % Ra
            for midx2 = 1
                nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 2nd-dimensional "num_nonzero34" of Ra_mn
                Ra(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn(:, idx),2);
            end
            % Rb
            if ~isempty(hbw)
                for midx2 = 1
                    nz_midx1 = find( SRb_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRb_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRb_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRb_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRb_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of Rb_mn
                    Rb(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn(:, idx),2);
                end
            end
        end
    end
    if sim.include_Raman
        clearvars Ra_mn;
        if ~isempty(hbw)
            clearvars Rb_mn;
        end
    end
end

% Calculate h*Ra as F-1(h F(Ra))
% The convolution using Fourier transform is faster if both arrays are
% large. If one of the array is small, "conv" can be faster.
% Please refer to
% "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
% for more information.
if sim.include_Raman
    Ra = fft(haw.*ifft(Ra,[],1),[],1);
    
    if ~isempty(hbw) % polarized fields with an anisotropic Raman
        Rb = fft(hbw.*ifft(Rb,[],1),[],1);
    end
    
    if isempty(hbw)
        nonlinear = Kerr + sum(Ra.*permute(At_wNoise,[1 3 2]),3);
    else % polarized fields with an anisotropic Raman
        nonlinear = Kerr + sum((Ra+Rb).*permute(At_wNoise,[1 3 2]),3);
    end
else
    nonlinear = Kerr;
end

% Gain term
func = solve_gain_rate_eqn_helpers();
gAt = func.solve_Power('signal',...
                       gain_rate_eqn,...
                       [],...
                       gain_rate_eqn.extended_cross_sections{2},...
                       N,...
                       [],[],Aw);

% Now everything has been summed into nonlinear, so transform into the
% frequency domain for the prefactor. It's transformed into the time domain
% for adding the gain term, following by tranforming back to the frequency
% domain for the output
dAwdz = ifft(fft(n2_prefactor.*ifft(nonlinear,[],1),[],1) + gAt,[],1);

end

%% RK4 for ASE
function A1t = RK4_ASE(A0w,...
                       gain_rate_eqn,...
                       dz,...
                       N,...
                       dAt_SE)

a1 = N_op_ASE(A0w,...
              gain_rate_eqn,...
              N);
a2 = N_op_ASE(A0w+a1*(dz/2),...
              gain_rate_eqn,...
              N);
a3 = N_op_ASE(A0w+a2*(dz/2),...
              gain_rate_eqn,...
              N);
a4 = N_op_ASE(A0w+a3*(dz),...
              gain_rate_eqn,...
              N);

A1 = A0w + (a1+2*a2+2*a3 +a4)*(dz/6);

A1t = fft(A1,[],1) + dAt_SE; % add spontaneous emission

end
% and the nonlinear operator in RK4_ASE
function dAwdz = N_op_ASE(Aw,...
                          gain_rate_eqn,...
                          N)
%N_op_ASE Calculate dAdz

% Gain term
func = solve_gain_rate_eqn_helpers();
gAt = func.solve_Power('signal',...
                       gain_rate_eqn,...
                       [],...
                       gain_rate_eqn.extended_cross_sections{1},...
                       N,...
                       [],[],Aw);

dAwdz = ifft(gAt,[],1);

end