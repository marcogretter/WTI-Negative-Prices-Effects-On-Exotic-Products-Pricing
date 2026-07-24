function otmSurface = buildFwdOtmSurfacePoint2(snapshot, futureExpiries, discountFactors, F_curve)
% The function returns otmSurface, a struct in which there are the
% prices of Fwd-OTM calls and puts.
%
% The output is a flat struct:
% each field is a column vector, and each row corresponds to one OTM option.

% Initialize the output
otmSurface.MaturityIndex = [];
otmSurface.Expiry = [];
otmSurface.TTM = [];
otmSurface.Strike = [];
otmSurface.Forward = [];
otmSurface.DiscountFactor = [];
otmSurface.Moneyness = [];
otmSurface.IsCall = [];
otmSurface.MarketPrice = [];
otmSurface.sigmaATM = [];
otmSurface.allCall = [];
otmSurface.chi = [];

valuationDate = snapshot(1).valuationDate;

for i = 1:length(futureExpiries)

    T = yearfrac(valuationDate, futureExpiries(i), 3);

    if T <= 0
        continue;
    end

    F = F_curve(i);
    B = discountFactors(i);

    if ~isfinite(F) || ~isfinite(B) || B <= 0
        continue;
    end

    
    % Calls OTM
    

    K_call = snapshot(i).K_call(:);
    C_call = snapshot(i).C_call(:);

    keepCall = K_call >= F;

    nCall = sum(keepCall);

    if nCall > 0
        otmSurface.MaturityIndex = [otmSurface.MaturityIndex; repmat(i, nCall, 1)];
        otmSurface.Expiry = [otmSurface.Expiry; repmat(futureExpiries(i), nCall, 1)];
        otmSurface.TTM = [otmSurface.TTM; repmat(T, nCall, 1)];
        otmSurface.Strike = [otmSurface.Strike; K_call(keepCall)];
        otmSurface.Forward = [otmSurface.Forward; repmat(F, nCall, 1)];
        otmSurface.DiscountFactor = [otmSurface.DiscountFactor; repmat(B, nCall, 1)];
        otmSurface.Moneyness = [otmSurface.Moneyness; K_call(keepCall) - F];
        otmSurface.IsCall = [otmSurface.IsCall; true(nCall, 1)];
        otmSurface.MarketPrice = [otmSurface.MarketPrice; C_call(keepCall)];
    end

    % Puts OTM
    
    K_put = snapshot(i).K_put(:);
    P_put = snapshot(i).P_put(:);

    keepPut = K_put < F;

    nPut = sum(keepPut);

    if nPut > 0
        otmSurface.MaturityIndex = [otmSurface.MaturityIndex; repmat(i, nPut, 1)];
        otmSurface.Expiry = [otmSurface.Expiry; repmat(futureExpiries(i), nPut, 1)];
        otmSurface.TTM = [otmSurface.TTM; repmat(T, nPut, 1)];
        otmSurface.Strike = [otmSurface.Strike; K_put(keepPut)];
        otmSurface.Forward = [otmSurface.Forward; repmat(F, nPut, 1)];
        otmSurface.DiscountFactor = [otmSurface.DiscountFactor; repmat(B, nPut, 1)];
        otmSurface.Moneyness = [otmSurface.Moneyness; K_put(keepPut) - F];
        otmSurface.IsCall = [otmSurface.IsCall; false(nPut, 1)];
        otmSurface.MarketPrice = [otmSurface.MarketPrice; P_put(keepPut)];
    end

end

end