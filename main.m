% Alessandro Di Nola
% University of Konstanz

clear
clc
close all

%% Set parameters

[par,vfi_options,sim_options] = set_params();

%% Make grids

par = make_grids(par);

%% Solve the model
if vfi_options.method == 1
    [sol,exit_flag] = sub_vfi(par,vfi_options);
elseif vfi_options.method == 2
    [sol,exit_flag] = sub_vfi_vec(par,vfi_options);
end

figure
plot(par.k_grid,sol.V1(:,1),'linewidth',2)
hold on
plot(par.k_grid,sol.V1(:,round(par.nz/2)),'linewidth',2)
hold on
plot(par.k_grid,sol.V1(:,end),'linewidth',2)
legend('zmin','zmed','zmax','location','best')
title("Value function")
xlabel("Capital today")
print('valfun','-dpng')

figure
plot(par.k_grid,par.k_grid,'--','linewidth',2)
hold on
plot(par.k_grid,sol.kpol(:,1),'linewidth',2)
hold on
plot(par.k_grid,sol.kpol(:,round(par.nz/2)),'linewidth',2)
hold on
plot(par.k_grid,sol.kpol(:,end),'linewidth',2)
legend('45 line','zmin','zmed','zmax','location','best')
title("Policy function")
xlabel("Capital today")
ylabel("Capital tomorrow")
print('polfun','-dpng')

%% Simulate panel data
sim = simmodel(par,sol,sim_options);

figure
plot(sim.k(24,:))
title('Simulated time series for capital')

figure
plot(sim.I(24,:))
title('Simulated time series for investment')

%% Compute model moments

[mom] = makemoments(par,sim);

% Save results
% if vfi_options.howard==0
%     V_nohoward = sol.V1;
%     save nohoward V_nohoward
% elseif vfi_options.howard==1
%     V_howard = sol.V1;
%     save howard V_howard
% end







