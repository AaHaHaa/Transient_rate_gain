function [x,fval,exitflag,output] = myfminsearch(funfcn,x,options,usual_delta,zero_term_delta)
%FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%
% -------------------------------------------------------------------------
% *Note:
%   This is copied and modified from MATLAB's (2024a) fminsearch internal function.
%   I removed redundant functions that aren't used in optimizing population differences.
%   These modifications are crucial to have the code work for different versions of MATLAB.
%
%   Most important of all,
%   I modified its convergence criteria to conform with the population optimization process. Check line 345-351.
%
%   MATLAB default stopping criteria of the Nelder-Mead method are
%   1. Values from all points are close to each other
%   2. Positions of all points are close to each other
%
%   Since I don't need all points to be close enough, which is a strict stopping criterion for our use,
%   I changed the stopping criteria into
%   1. (Same as default) Values from all points are close to each other
%   2. The minimum value from the first point of the Nelder-Mead method is good enough (smaller than the preset threshold)
%      This best point is chosen as the output that is able to give the minimum value for the optimization problem.
%
%   Besides, two inputs are added:
%       usual_delta: 5 percent deltas for non-zero terms
%       zero_term_delta: Even smaller delta for zero elements of x
%
%       by Yi-Hao Chen, PhD in Frank Wise's group at Applied Physics, Cornell University, 1/17/2025
% -------------------------------------------------------------------------
%
%   Below are retained from MATLAB's original comments:
%
%   X = FMINSEARCH(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X and 
%   returns a scalar function value F evaluated at X. X0 can be a scalar, vector 
%   or matrix.
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, FunValCheck,
%   PlotFcns, and OutputFcn.
%
%   X = FMINSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminsearch' in PROBLEM.solver. 
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns an EXITFLAG that describes
%   the exit condition. Possible values of EXITFLAG and the corresponding
%   exit conditions are
%
%    1  Maximum coordinate difference between current best point and other
%       points in simplex is less than or equal to TolX, and corresponding 
%       difference in function values is less than or equal to TolFun.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value
%     SIN evaluated at X.
%
%     FUN can be an anonymous function:
%        X = fminsearch(@(x) norm(x),[1;2;3])
%     returns a point near the minimizer [0;0;0].
%
%     FUN can be a parameterized function. Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
%        c = 1.5;                        % The parameter.
%        X = fminsearch(@(x) f(x,c),[0.3;1])
%        
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

%   Copyright 1984-2023 The MathWorks, Inc.

n = numel(x);

%% Set the optimization parameters
defaultopt = makeDefaultopt();
optimgetFlag = 'fast';
tolx = optimget(options,'TolX',defaultopt,optimgetFlag); %#ok
tolf = optimget(options,'TolFun',defaultopt,optimgetFlag);
maxfun = optimget(options,'MaxFunEvals',defaultopt,optimgetFlag);
maxiter = optimget(options,'MaxIter',defaultopt,optimgetFlag);

if ischar(maxfun) || isstring(maxfun)
    if strcmpi(maxfun,'200*numberofvariables')
        maxfun = 200*n;
    else
        error('MATLAB:fminsearch:OptMaxFunEvalsNotInteger',...
            getString(message('MATLAB:optimfun:fminsearch:OptMaxFunEvalsNotInteger')));
    end
end
if ischar(maxiter) || isstring(maxiter)
    if strcmpi(maxiter,'200*numberofvariables')
        maxiter = 200*n;
    else
        error('MATLAB:fminsearch:OptMaxIterNotInteger',...
            getString(message('MATLAB:optimfun:fminsearch:OptMaxIterNotInteger')));
    end
end

%% Initialize parameters
rho = 1; 
chi = 2; 
psi = 0.5; 
sigma = 0.5;
np1 = n + 1;

%% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,np1); 
fv = zeros(1,np1);
v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = xin;    % Change x to the form expected by funfcn
fv(:,1) = funfcn(x);
itercount = 0;
% Initial simplex setup continues later

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
%usual_delta = 0.05;             % 5 percent deltas for non-zero terms
%zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
for j = 1:n
    y = xin;
    if y(j) ~= 0
        y(j) = (1 + usual_delta)*y(j);
    else
        y(j) = zero_term_delta;
    end
    v(:,j+1) = y;
    x(:) = y; 
    f = funfcn(x);
    fv(1,j+1) = f;
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv);
v = v(:,j);

itercount = itercount + 1;
func_evals = np1;

%% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v3-v1||,...,||v(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations 
% are exceeded
while func_evals < maxfun && itercount < maxiter
    %if all(abs(fv(1)-fv(2:np1)) <= max(tolf,10*eps(fv(1)))) && ...
    %        all(abs(v(:,2:np1)-v(:,1)) <= max(tolx,10*eps(max(v(:,1)))),'all')
    %    break
    %end
    if all(abs(fv(1)-fv(2:np1)) <= max(tolf/1e2,10*eps(fv(1)))) || fv(1) < tolf % stop if the best point is good enough
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the n (NOT n+1) best points
    xbar = sum(v(:,1:n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,np1);
    x(:) = xr; 
    fxr = funfcn(x);
    func_evals = func_evals+1;
    
    if fxr < fv(1)
        % Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,np1);
        x(:) = xe; 
        fxe = funfcn(x);
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,np1) = xe;
            fv(np1) = fxe;
            how = 'expand'; %#ok
        else
            v(:,np1) = xr;
            fv(np1) = fxr;
            how = 'reflect'; %#ok
        end
    else % fv(:,1) <= fxr
        if fxr < fv(n)
            v(:,np1) = xr;
            fv(np1) = fxr;
            how = 'reflect'; %#ok
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(np1)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,np1);
                x(:) = xc; 
                fxc = funfcn(x);
                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,np1) = xc;
                    fv(np1) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,np1);
                x(:) = xcc; 
                fxcc = funfcn(x);
                func_evals = func_evals+1;
                
                if fxcc < fv(np1)
                    v(:,np1) = xcc;
                    fv(np1) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j = 2:np1
                    v(:,j) = v(:,1)+sigma*(v(:,j) - v(:,1));
                    x(:)   = v(:,j); 
                    fv(j)  = funfcn(x);
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fv,j] = sort(fv);
    v = v(:,j);
    itercount = itercount + 1;
end   % while

x(:) = v(:,1);
fval = fv(:,1);

if func_evals >= maxfun
    exitflag = 0;
elseif itercount >= maxiter
    exitflag = 0;
else
    exitflag = 1;
end

output.iterations = itercount;
output.funcCount = func_evals;
output.algorithm = 'Nelder-Mead simplex direct search';

%--------------------------------------------------------------------------
function defaultopt = makeDefaultopt()
defaultopt = struct('MaxIter','200*numberOfVariables',...
                    'MaxFunEvals','200*numberOfVariables',...
                    'TolX',1e-4,...
                    'TolFun',1e-4);
