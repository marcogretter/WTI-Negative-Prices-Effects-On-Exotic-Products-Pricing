function [sigmaATMByMaturity,otmSurface] = estimateSigmaATMByMaturity(otmSurface)
% This function implements the second part of calibration of our project in
% which from the price of Calls and Puts OTM we have to calibrate the sigma
% ATM

% We search for the maximum index in the maturities (9 in our case)
maxMatIdx = max(otmSurface.MaturityIndex);
% Initialize vectors of sigma_ATM
sigmaATMByMaturity = nan(maxMatIdx, 1);

nRows = numel(otmSurface.MarketPrice);

otmSurface.allCall = nan(nRows, 1);
otmSurface.sigmaATM = nan(nRows, 1);

% We get the indices of the maturities that we have and we mantain the
% original order
uniqueMaturities = unique(otmSurface.MaturityIndex, 'stable');

for j = 1:numel(uniqueMaturities)
    
    matIdx = uniqueMaturities(j);

    % We get only the lines related to that maturity index
    rows = otmSurface.MaturityIndex == matIdx;

    x = otmSurface.Moneyness(rows);
    B = otmSurface.DiscountFactor(rows);
    T = otmSurface.TTM(rows);

    % We see for every row if we are handling a call
    isCall = otmSurface.IsCall(rows);
    marketPrice = otmSurface.MarketPrice(rows);
    
    callEqPrice = marketPrice;
    % By put-call parity we know: C-P = -Bx
    % IMP=> marketPrice(~isCall) is the price of a put, since we have the
    % formula from the paper only for the put, we have to reconduct
    % ourselves to the price of a call
    callEqPrice(~isCall) = marketPrice(~isCall) - B(~isCall) .* x(~isCall);
    otmSurface.allCall(rows) = callEqPrice;
    
    % We get how many options we have for the j-th maturity
    n = numel(x);
    impVol = nan(n, 1);

    % You get the Bachelier implied volatility
    for i = 1:n
        impVol(i) = bachelierImpliedVolFromCall(x(i), B(i), T(i), callEqPrice(i));
    end

    % We mantain only the valid impl_volatilities
    valid = isfinite(impVol) & impVol > 0;

    xValid = x(valid);
    volValid = impVol(valid);

    % Check if we have at least a valid imp_vol
    if numel(volValid) < 1
        continue;
    end

    % We order the imp_vols putting as first the closest one to the ATM (x=0) 
    [~, ord] = sort(abs(xValid), 'ascend');
    % We take the 7 implied_vols closest to x=0
    nNear = min(7, numel(ord));
    % We get the related indices
    sel = ord(1:nNear);

    
    xNear = xValid(sel);
    volNear = volValid(sel);

    % If I have at least 3 impl_vols and at least one on both sides I do 
    % QUADRATIC interpolation to find the sigma_ATM, else I take the 
    % closest one (we take into account the smile)
    hasLeftWing = any(xNear < 0);
    hasRightWing = any(xNear > 0);
    
    if numel(volNear) >= 3 && hasLeftWing && hasRightWing
        p = polyfit(xNear, volNear, 2);
        sigmaATM = polyval(p, 0);
    else
        sigmaATM = volNear(1);
    end
    
    % We take only values of volatility > 0
    sigmaATMByMaturity(matIdx) = max(sigmaATM, 1e-8);
    otmSurface.sigmaATM(rows) = sigmaATMByMaturity(matIdx);
    otmSurface.sigmaATM = otmSurface.sigmaATM(:);
    otmSurface.chi(rows) = otmSurface.Moneyness(rows) ./ (otmSurface.sigmaATM(rows) .* sqrt(otmSurface.TTM(rows)));
    otmSurface.chi = otmSurface.chi(:);
end

end