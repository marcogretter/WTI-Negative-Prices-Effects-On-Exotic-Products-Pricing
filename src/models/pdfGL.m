function y = pdfGL(x, alpha, beta)
%PDFGL Generalized Logistic Type IV pdf with zero mean.
% alpha: left tail parameter
% beta: right tail parameter

gamma = psi(beta) - psi(alpha);
z = x - gamma;
logC = gammaln(alpha + beta) - gammaln(alpha) - gammaln(beta);
logY = logC + alpha .* z - (alpha + beta) .* log1p(exp(z));
y = exp(logY);

end