function [eta,kappa] = calibrate_eta_kappa(otmSurface, eta_0,kappa_0,M,value,type)

p0 = [eta_0, kappa_0];
options = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-6,'TolFun',1e-8);
p_results = fminsearch(@(p) obj_fun_AB(otmSurface,p(1),p(2),M,value,type), p0, options);
eta = p_results(1);
kappa = p_results(2);
end

