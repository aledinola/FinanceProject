function [par,vfi_options,sim_options] = set_params()

% Economic parameters
par.alpha  = 0.5;   % returns to scale production function
par.delta  = 0.06;  % depreciation rate
par.beta   = 0.96;  % discount factor
par.lambda = 0.05;  % linear cost of external finance
par.rho    = 0.75;  % autocorrelation of productivity shock
par.sigma  = 0.30;  % standard deviation of productivity shock

% Numerical parameters

% Grids
par.nz = 21;     % no. of grid points for productivity shock
par.nk = 1000;  % no. of grid points for capital

par.k_min = 1e-6;
par.k_max = 200;

% VFI properties
vfi_options.method  = 1;    % 1=loops, 2=vectorized
vfi_options.howard  = 1;    % flag 0-1
vfi_options.n_howard = 50;  % no. of iterations for Howard policy improvement
vfi_options.tol     = 1e-6; % tolerance
vfi_options.maxiter = 1000; % max no. of iterations
vfi_options.verbose = 1;    % 0-1-2 flag

% Simulation
sim_options.Tburn   = 1000;
sim_options.Tsim    = 100;
sim_options.Ttot    = sim_options.Tburn+sim_options.Tsim;
sim_options.Nfirms = 10000;
sim_options.verbose = 1;

end % end function <set_params>