function [sim] = simmodel(par,sol,sim_options)

% DESCRIPTION:
%	Simulate a panel dataset from the model
% INPUTS:
% 	"par" Structure with model parameters and grids
% 	"sol" Structure with model policy function
%
% OUTPUTS:
%	"sim" N*T simulated data

disp("Doing simulation...")

%% Unpacking:
Tsim  = sim_options.Tsim;
Tburn = sim_options.Tburn;
Ttot  = Tsim+Tburn;
N     = sim_options.Nfirms;

%kpol     = sol.kpol;
kpol_ind = sol.kpol_ind;

sim.k         = zeros(N,Ttot+1);
sim.k_ind     = ones(N,Ttot+1);
sim.y         = zeros(N,Ttot+1); % Revenues or profits
sim.I         = zeros(N,Ttot+1); % Investment
sim.cash_flow = zeros(N,Ttot+1); % Cash flow (profit minus investment)
sim.Irate     = zeros(N,Ttot+1);
sim.prof      = zeros(N,Ttot+1); % Profitability: revenues/capital

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

sim.k_ind(:,1) = round(par.nk/2); % Initial condition for capital
% for i = 1:N
%     for t = 1:Ttot+1
%         if t<=Ttot
%             % Use policy function to simulate endogenous state (indexes)
%             x = sim.k_ind(i,t)+(z_sim_ind(i,t)-1)*par.nk;
%             sim.k_ind(i,t+1) = kpol_ind(x);
%             %sim.k_ind(i,t+1) = kpol_ind(sim.k_ind(i,t),z_sim_ind(i,t));
%         end
%         sim.k(i,t) = par.k_grid(sim.k_ind(i,t)); % capital in values
%     end
% end

for t = 1:Ttot+1
    if t<=Ttot
        % Use policy function to simulate endogenous state (indexes)
        % sub2ind()
        sim.k_ind(:,t+1) = kpol_ind(sim.k_ind(:,t)+(z_sim_ind(:,t)-1)*par.nk);
    end
    sim.k(:,t) = par.k_grid(sim.k_ind(:,t)); % capital in values
end

%% Simulate other variables

%for i = 1:N
    for t = 1:Ttot+1
        if t<=Ttot
            % investment
            sim.I(:,t+1) = sim.k(:,t+1)-(1-par.delta)*sim.k(:,t);
            % investment rate is (k'-(1-delta)k)/k
            sim.Irate(:,t+1) = sim.I(:,t+1)./sim.k(:,t);
        end
        % output
        sim.y(:,t) = fun_prod(sim.k(:,t),z_sim(:,t),par);
        % cash flow
        sim.cash_flow(:,t) = sim.y(:,t)-sim.I(:,t);
        %sim.Irate(:,t) = sim.I(:,t)./sim.k(:,t); % investment rate
        sim.prof(:,t) = sim.y(:,t)./sim.k(:,t);
        % add external finance or equity issuances
    end
%end

%% Discard burn-in

keep = (Tburn+1):(Ttot-1);
sim.k = sim.k(:,keep);
sim.I = sim.I(:,keep);
sim.y = sim.y(:,keep);
sim.cash_flow = sim.cash_flow(:,keep);
sim.Irate = sim.Irate(:,keep);

if sim_options.verbose>=1
    disp("Runtime for simulation:")
    toc;
end

%% Simulation checks

fprintf(' \n')
disp("Simulation checks:")
fprintf('Minimum value of ksim = %f \n', min(sim.k(:)))
fprintf('Maximum value of ksim = %f \n', max(sim.k(:)))
fprintf('In-sample average of ksim = %f \n', mean(sim.k(:)))

end %end function <simmodel>



