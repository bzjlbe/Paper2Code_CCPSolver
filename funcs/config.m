function params = config(rho, h, num, e_imp, tol)
% Parameter configuration file
% Input parameters:
%   rho: generalized_alpha method parameter
%   h: time step size
%   num: number of steps
%   e_imp: impact restitution coefficient
%   tol: tolerance

% Generalized-alpha method parameters
params.alpha_m = (2*rho-1)/(1+rho);
params.alpha_f = rho/(1+rho);
params.gama = 0.5-params.alpha_m+params.alpha_f;
params.eta = (0.5+params.gama)^2/4;
params.eta_t = h*params.eta/params.gama;
params.gama_t = (1-params.alpha_m)/(h*params.gama*(1-params.alpha_f));

% Simulation parameters
params.h = h;
params.num = num;
params.time = (0:h:h*num)';
params.e_imp = e_imp;
params.tol = tol;
end