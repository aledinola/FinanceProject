function [sol,exit_flag] = sub_vfi_vec(par,vfi_options)

% DESCRIPTION:
%	Solve the model with discrete value function iteration
%   Wrt <sub_vfi>, I vectorize the loop over k
% INPUTS:
% 	"par" Structure with model parameters and grids
%
% OUTPUTS:
%	"sol" Structure with value function V(k,z) and policy function k'(k,z)
%   "exit_flag" Flag to indicate convergence (a negative value means lack of convergence)

% Unpack some useful variables
nk     = par.nk;
nz     = par.nz;
beta   = par.beta;
k_grid = par.k_grid;
z_grid = par.z_grid;
z_prob = par.z_prob;

disp("Start VFI, vectorized")

% Precompute current-period payoff E(k',k,z)
if vfi_options.verbose>=1
    tic;
end
payoff = zeros(nk,nk,nz);

for z_c = 1:nz
    for k_c = 1:nk
        k_val  = k_grid(k_c);
        z_val  = z_grid(z_c);
        kp_val = k_grid;
        payoff(:,k_c,z_c) = ReturnFn(kp_val,k_val,z_val,par);
    end
end

if vfi_options.verbose>=1
    disp("Time to precompute static payoff")
    toc;
    fprintf(" \n")
end

% Initialize VFI

iter = 1;
err  = 1+vfi_options.tol;
V0   = zeros(nk,nz);
V1   = zeros(nk,nz);
kpol = zeros(nk,nz);
exit_flag = 0;
if vfi_options.verbose>=1
    tic;
end

while (err>vfi_options.tol && iter<=vfi_options.maxiter)
    
    % Compute expected value function EV(k',z)
    EV = V0*z_prob'; % dim: (k',z)
    
    % Maximization step to generate V1 - vectorized
    % RHS_vec is a matrix with (k',k) i.e. choices on the 1dim,states on
    % the 2nd dim
    for z_c = 1:nz
        RHS_vec           = payoff(:,:,z_c)+beta*EV(:,z_c);
        [V1(:,z_c),k_ind] = max(RHS_vec);
        kpol(:,z_c)       = k_grid(k_ind);
    end
    
    err = max(abs(V1(:)-V0(:)));
    
    if vfi_options.verbose>=1
        fprintf("iter = %d, err = %f \n",iter,err)
    end
    
    iter = iter+1;
    V0 = V1;
    
end % end while

if vfi_options.verbose>=1
    disp("Runtime for VFI:")
    toc;
end

if err>vfi_options.tol
    warning("VFI did not converge!")
    exit_flag=-1;
end

% Pack results

sol.V1 = V1;
sol.kpol = kpol;

end % end function <sub_vfi_vec>





