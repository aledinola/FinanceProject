function [sol,exit_flag] = sub_vfi(par,vfi_options)

% DESCRIPTION:
%	Solve the model with discrete value function iteration
%   Check out also "sub_vfi_vec" for a faster(?) implementation
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
n_howard = vfi_options.n_howard;

disp("Start VFI with loops")

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
kpol_ind = ones(nk,nz);
exit_flag = 0;
if vfi_options.verbose>=1
    tic;
end

while (err>vfi_options.tol && iter<=vfi_options.maxiter)
    
    % Compute expected value function EV(k',z)
    EV = V0*z_prob'; % dim: (k',z)
    
    % Maximization step to generate V1 - with loops
    for z_c = 1:nz
        for k_c = 1:nk
            RHS_vec             = payoff(:,k_c,z_c)+beta*EV(:,z_c);
            [V1(k_c,z_c),k_ind] = max(RHS_vec);
            kpol_ind(k_c,z_c)   = k_ind;
            kpol(k_c,z_c)       = k_grid(k_ind);
        end
    end
    
    %-------------- Howard ----------------------------%
    if vfi_options.howard==1
        payoff_max = zeros(nk,nk);
        for z_c = 1:nz
            for k_c = 1:nk
                payoff_max(k_c,z_c) = payoff(kpol_ind(k_c,z_c),k_c,z_c);
            end
        end
        
        for h_c = 1:n_howard % Howard iteration counter
            for z_c = 1:nz
                %for k_c = 1:nk
                    %EVh = dot(z_prob(z_c,:),V1(kpol_ind(k_c,z_c),:));
                    EVh = V1(kpol_ind(:,z_c),:)*z_prob(z_c,:)';
                    V1(:,z_c) = payoff_max(:,z_c)+beta*EVh;
                %end
            end
        end
        
    end
    %-------------- End Howard ----------------------------%
    
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

disp("VFI diagnostics: make sure grid is not too restrictive")
fprintf('Lower bound capital grid = %f \n',par.k_grid(1))
fprintf('Min of capital policy    = %f \n',min(kpol(:)))
fprintf('Upper bound capital grid = %f \n',par.k_grid(end))
fprintf('Max of capital policy    = %f \n',max(kpol(:)))
fprintf(' \n')

% Pack results

sol.V1 = V1;
sol.kpol = kpol;
sol.kpol_ind = kpol_ind;

end % end function <sub_vfi>
