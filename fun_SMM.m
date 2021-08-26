function [obj_smm,mom] = fun_SMM(par,vfi_options,sim_options)


%% Make grids

par = make_grids(par);

%% Solve the model
if vfi_options.method == 1
    [sol,exit_flag] = sub_vfi(par,vfi_options);
elseif vfi_options.method == 2
    [sol,exit_flag] = sub_vfi_vec(par,vfi_options);
end

if exit_flag<0
    warning('VFI did not converge')
    obj_smm = +inf; mom = [];
    return
end

%% Simulate panel data
sim = simmodel(par,sol,sim_options);

%% Compute model moments

[mom] = makemoments(par,sim);

obj_smm = [];

end % end function <fun_SMM>

