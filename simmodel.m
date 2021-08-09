function [sim] = simmodel(par,sol,sim_options)

% DESCRIPTION:
%	Simulate a panel dataset from the model
% INPUTS:
% 	"par" Structure with model parameters and grids
% 	"sol" Structure with model policy function
%
% OUTPUTS:
%	"sim" N*T simulated data

%% Unpacking:
Tsim = sim_options.Tsim;
Tburn = sim_options.Tburn;
Ttot = Tsim+Tburn;
N = sim_options.Nfirms;

kpol     = sol.kpol;
kpol_ind = sol.kpol_ind;

sim.k         = zeros(N,Ttot+1);
sim.k_ind     = ones(N,Ttot+1);
sim.y         = zeros(N,Ttot+1); % Revenues or profits
sim.I         = zeros(N,Ttot+1); % Investment
sim.cash_flow = zeros(N,Ttot+1); % Cash flow (profit minus investment)

if sim_options.verbose>=1
    tic;
end

%% Draw uniform shocks in (0,1)
u = rand(N,Ttot+1); % N*Ttot matrix of (0,1) random numbers

%% Obtain z_sim_ind

% Given u and given the Markov chain for the exogenous stochastic variables
% z (productivity shock), simulate a panel of realizations for z

% idxM = markov_sim(nInd, T, prob0V, trProbM, rvInM, dbg)
prob0V    = ones(par.nz,1)/par.nz;
z_sim_ind = markov_sim(N, Ttot+1, prob0V, par.z_prob', u, 1);
z_sim = par.z_grid(z_sim_ind);

%% Simulate endogenous state

sim.k_ind(:,1) = round(par.nk/2);
for i = 1:N
    for t = 1:Ttot+1
        if t<=Ttot
            % Use policy function to simulate endogenous state (indexes)
            sim.k_ind(i,t+1) = kpol_ind(sim.k_ind(i,t),z_sim_ind(i,t));
        end
        sim.k(i,t) = par.k_grid(sim.k_ind(i,t)); % capital in values
        
    end
end

%% Simulate other variables

for i = 1:N
    for t = 1:Ttot+1
        sim.y(i,t) = fun_prod(sim.k(i,t),z_sim(i,t),par);
        if t<=Ttot
            sim.I(i,t+1) = sim.k(i,t+1)-(1-par.delta)*sim.k(i,t);
        end
        sim.cash_flow(i,t) = sim.y(i,t)-sim.I(i,t);
    end
end

%% Discard burn-in

keep = (Tburn+1):(Ttot-1);
sim.k = sim.k(:,keep);
sim.I = sim.I(:,keep);
sim.y = sim.y(:,keep);
sim.cash_flow = sim.cash_flow(:,keep);

if sim_options.verbose>=1
    disp("Runtime for simulation:")
    toc;
end

%% Simulation checks

fprintf(' \n')
fprintf('Minimum value of ksim = %f \n', min(sim.k(:)))
fprintf('Maximum value of ksim = %f \n', max(sim.k(:)))

end %end function <simmodel>



