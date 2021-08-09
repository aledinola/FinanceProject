function F = fun_prod(k,z,par)

% DESCRIPTION:
%	Production function f(k,z) = z*k^alpha
% INPUTS:
%   "k"   capital
%   "z"   productivity shock
%   "par" structure with model parameters

F = z*k.^par.alpha;

end % end function <fun_prod>