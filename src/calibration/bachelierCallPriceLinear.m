function callPrice = bachelierCallPriceLinear(moneyness, discountFactor, normalVol, timeToMaturity)
% Closed-form Bachelier call in linear moneyness.
% equations 5,6 in the [3] Additive Bachelier paper.
%
% Inputs:
%   moneyness      : n x 1 vector, x = K - F(0,T)
%   discountFactor : n x 1 vector, B(0,T)
%   normalVol      : n x 1 vector, Bachelier normal volatility
%   timeToMaturity : n x 1 vector, maturity in years
%
% Output:
%   callPrice      : n x 1 vector

x = moneyness(:);

stdDev = normalVol(:) .* sqrt(timeToMaturity(:));

y = x ./ stdDev;

c_b = @(y, sigma) -y .* normcdf(-y ./ sigma) + sigma .* normpdf(-y ./ sigma);

callPrice = discountFactor(:) .* stdDev .* c_b(y, 1);

end



