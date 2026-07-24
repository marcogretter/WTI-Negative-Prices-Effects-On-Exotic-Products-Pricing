function impliedVol = bachelierImpliedVolFromCall(x, B, T, callPrice)

intrinsic = B * max(-x, 0);

if callPrice <= intrinsic
    impliedVol = NaN;
    return;
end

priceError = @(vol) bachelierCallPriceLinear(x, B, vol, T) - callPrice;

volLow = 1e-8;
volHigh = 1;

while priceError(volHigh) < 0
    volHigh = 2 * volHigh;

    if volHigh > 1e4
        impliedVol = NaN;
        return;
    end
end

impliedVol = fzero(priceError, [volLow, volHigh]);

end