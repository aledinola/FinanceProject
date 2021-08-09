function F = ReturnFn(kp_val,k_val,z_val,par)

% DESCRIPTION:
%	Current-period payoff E(k',k,z) to shareholders
% INPUTS:
%   "kp_val" k',next-period capital
%   "k_val"  current period capital
%   "z_val"  current period productivity shock
% OUTPUTS:
%   "F"      Distribution to shareholders, if negative is equity injection

investment = kp_val-(1-par.delta)*k_val; % I=k'-(1-delta)k
cash_flow  = fun_prod(k_val,z_val,par)-investment;

% if cash_flow>=0
%     F = cash_flow;
% else
%     F = (1+par.lambda)*cash_flow;
% end
F = cash_flow;
F(cash_flow<0) = F(cash_flow<0)+ par.lambda*cash_flow(cash_flow<0);

end % end function <ReturnFn>