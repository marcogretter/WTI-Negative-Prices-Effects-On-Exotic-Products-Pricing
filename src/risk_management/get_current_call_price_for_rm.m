function [currentCallPrice, matchedStrike, isExactMatch] = get_current_call_price_for_rm(snapshotMaturity, originalStrike)

callStrikes = snapshotMaturity.K_call(:);
callPrices  = snapshotMaturity.C_call(:);

if isempty(callStrikes) || isempty(callPrices)
    error('No call quotes available for the selected maturity.');
end

idxExact = find(abs(callStrikes - originalStrike) < 1e-12, 1);

if ~isempty(idxExact)
    matchedStrike = callStrikes(idxExact);
    currentCallPrice = callPrices(idxExact);
    isExactMatch = true;
else
    [~, idxNearest] = min(abs(callStrikes - originalStrike));
    matchedStrike = callStrikes(idxNearest);
    currentCallPrice = callPrices(idxNearest);
    isExactMatch = false;
end

end