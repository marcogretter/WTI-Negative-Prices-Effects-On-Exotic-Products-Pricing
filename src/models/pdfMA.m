function y = pdfMA(x, alpha, beta)
% PDFMA Asymmetric Laplace pdf with zero mean.
%

coeff = (alpha * beta)/ (alpha + beta);
gamma = 1 / alpha - 1 / beta;

left = x < gamma;
right = ~left;

y(left) = coeff * exp(alpha .* (x(left) - gamma));
y(right) = coeff * exp(-beta .* (x(right) - gamma));

end 