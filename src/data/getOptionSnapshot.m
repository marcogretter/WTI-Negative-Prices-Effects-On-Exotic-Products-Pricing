function snapshot = getOptionSnapshot(data, valuationDate)
%GETOPTIONSNAPSHOT Extract option data for one valuation date.
%
% Output for each maturity:
%
%   Point 1:
%       K, C, P
%       only common strikes with valid call and put prices.
%
%   Point 2:
%       K_call, C_call
%       all valid calls.
%
%       K_put, P_put
%       all valid puts.

snapshot = [];

for j = 1:numel(data)
    row = find(data(j).dates == valuationDate, 1);

    if isempty(row)
        warning('Valuation date not found for maturity %s. Skipping.', ...
        data(j).maturityCode);
        continue;
    end

    % Raw call data
    
    K_call_all = data(j).callStrikes(:);
    C_all = data(j).calls(row, :).';

    
    % Raw put data
    K_put_all = data(j).putStrikes(:);
    P_all = data(j).puts(row, :).';

    
    % Point 2 data: keep calls and puts separately
    
    validCall = isfinite(K_call_all) & isfinite(C_all);
    validPut  = isfinite(K_put_all)  & isfinite(P_all);

    K_call = K_call_all(validCall);
    C_call = C_all(validCall);

    K_put = K_put_all(validPut);
    P_put = P_all(validPut);

    
    % Point 1 data: common strikes only
    
    [K_common_raw, idxCall, idxPut] = intersect(K_call_all, K_put_all, 'stable');

    C_common_raw = C_all(idxCall);
    P_common_raw = P_all(idxPut);

    validCommon = isfinite(K_common_raw) & ...
                      isfinite(C_common_raw) & ...
                      isfinite(P_common_raw);

    K_common = K_common_raw(validCommon);
    C_common = C_common_raw(validCommon);
    P_common = P_common_raw(validCommon);
    % Save in snapshot
    
    item = struct();

    item.maturityCode = data(j).maturityCode;
    item.valuationDate = valuationDate;

    % Data for Point 1: put-call parity / synthetic forward
    item.K = K_common;
    item.C = C_common;
    item.P = P_common;

    % Data for Point 2: Fwd-OTM surface
    item.K_call = K_call;
    item.C_call = C_call;

    item.K_put = K_put;
    item.P_put = P_put;

    if isempty(snapshot)
        snapshot = item;
    else
        snapshot(end+1) = item; %
    end

end

end