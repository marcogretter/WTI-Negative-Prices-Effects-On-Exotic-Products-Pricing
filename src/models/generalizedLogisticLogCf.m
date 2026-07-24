function logPhi = generalizedLogisticLogCf(u, alpha, betaPar)

gammaShift = psi(betaPar) - psi(alpha);

logPhi = 1i .* u .* gammaShift + loggamma(alpha + 1i .* u) ...
    + loggamma(betaPar - 1i .* u) - gammaln(alpha) - gammaln(betaPar);

end

function [f] = loggamma(z)
% LOGGAMMA  Log-Gamma function valid in the complex plane.
%           Lanczos approximation in logarithmic form.
%
%           z may be complex and of any size.
%
% Usage:
%   f = loggamma(z)
%
% Output:
%   f = log(Gamma(z))

siz = size(z);
z = z(:);
zz = z;

f = 0 .* z; % reserve space in advance

p = find(real(z) < 0);
if ~isempty(p)
   z(p) = -z(p);
end

g = 607/128;

c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];

z = z - 1;
zh = z + 0.5;
zgh = zh + g;

ss = 0.0;
for pp = size(c,1)-1:-1:1
    ss = ss + c(pp+1) ./ (z + pp);
end

sq2pi = 2.5066282746310005024157652848110;

% Log-Lanczos formula:
% Gamma(z+1) = sqrt(2*pi) * (c1 + ss) * (z+g+1/2)^(z+1/2) * exp(-(z+g+1/2))
f = log(sq2pi) + log(c(1) + ss) + zh .* log(zgh) - zgh;

% Gamma(1) = Gamma(2) = 1, therefore log Gamma = 0
f(z == 0 | z == 1) = 0.0;

% Adjust for negative real parts using reflection formula
if ~isempty(p)
   f(p) = log(-pi ./ (zz(p) .* sin(pi .* zz(p)))) - f(p);
end

% Adjust for negative poles
poles = find(round(zz) == zz & imag(zz) == 0 & real(zz) <= 0);
if ~isempty(poles)
   f(poles) = Inf;
end

f = reshape(f, siz);

end