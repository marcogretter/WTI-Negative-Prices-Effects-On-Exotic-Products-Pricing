function phi = cfAB(u, eta, k, alpha, scale)
%CFAB Characteristic function of the Additive Bachelier model.
%
% Inputs:
%   u     : Fourier variable
%   eta   : skew parameter
%   k     : tail / vol-of-vol parameter
%   alpha : stability parameter, alpha = 1/2 in the assignment
%   scale : sigma_t * sqrt(t)
%
% Model:
%   log phi(u) = psi(i*u*eta*scale + 0.5*u^2*scale^2; k, alpha)
%                + i*u*eta*scale

    v = 1i .* u .* eta .* scale + 0.5 .* (u.^2) .* scale.^2;

    if alpha == 0
        levyExponent = -(1 ./ k) .* log(1 + k .* v);
    else
        levyExponent = ((1 - alpha) ./ (alpha .* k)) .* ...
            (1 - (1 + k .* v ./ (1 - alpha)).^alpha);
    end

    phi = exp(levyExponent + 1i .* u .* eta .* scale);

end