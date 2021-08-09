function par = make_grids(par)

% Discretize productivity shock z via Tauchen's method
%tauchen(meanz, stdinnov, rho, multiple, znum)

[z_grid,par.z_prob] = tauchen(0,par.sigma,par.rho,3,par.nz);
par.z_grid = exp(z_grid');

% Grid for capital k
par.k_grid = linspace(par.k_min,par.k_max,par.nk)';


end % end function <make_grids>