function [price_exotic, diagnostics] = price_exotic_no_grid_from_increments( ...
    exoticType, x_T0_T1, x_T1_T2, ...
    B_T1, B_T2, F_T0_T2, K2, K1, num_grid_points)

%PRICE_EXOTIC_NO_GRID_FROM_INCREMENTS Numerical pricing
% for CoC, PoP and Chooser from already simulated additive increments.
%
% This version avoids both:
%   1) the stochastic mesh grid on X_T1;
%   2) the interpolation step.
%
% Assumption:
%   The simulated process is linear:
%       X_T2 = X_T1 + DeltaX_T1_T2
%
%   The terminal underlying is:
%       S_T2 = F(T0,T2) + X_T2
%
% Therefore:
%       S_T2 - K2 = X_T2 + F(T0,T2) - K2
%
% Inputs:
%   exoticType      : 'CoC', 'PoP', or 'Chooser'
%   x_T0_T1         : N_outer x 1 simulated increment from T0 to T1
%   x_T1_T2         : N_inner x 1 simulated increment from T1 to T2
%   B_T1, B_T2      : discount factors B(0,T1), B(0,T2)
%   F_T0_T2         : forward F(T0,T2)
%   K2              : strike of the inner vanilla option
%   K1              : strike of the compound option, ignored for Chooser
%   num_grid_points : unused, kept only for compatibility with mesh version
%
% Outputs:
%   price_exotic    : numerical exotic price
%   diagnostics     : struct with useful checks

if nargin < 9
    num_grid_points = [];
end

if nargin < 8 || isempty(K1)
    K1 = 0;
end

exoticType = validatestring(exoticType, {'CoC', 'PoP', 'Chooser'});

x_T1 = x_T0_T1(:);
dx_T1_T2 = x_T1_T2(:);

if isempty(x_T1) || isempty(dx_T1_T2)
    error('Input increment vectors must be non-empty.');
end

if any(~isfinite(x_T1)) || any(~isfinite(dx_T1_T2))
    error('Input increment vectors must contain only finite values.');
end

if ~isscalar(B_T1) || ~isscalar(B_T2) || B_T1 <= 0 || B_T2 <= 0
    error('Discount factors must be positive scalars.');
end

if ~isscalar(F_T0_T2) || ~isscalar(K2) || ~isscalar(K1)
    error('F_T0_T2, K2 and K1 must be scalar values.');
end

B12 = B_T2 / B_T1;

% Moneyness adjustment:
% S_T2 - K2 = X_T2 + F_T0_T2 - K2
m2 = F_T0_T2 - K2;

%% Inner vanilla values without grid/interpolation
%
% For each outer scenario x_T1(i), we need:
%
%   call_in_T1(i) = E[(x_T1(i) + DeltaX_T1_T2 + m2)^+]
%   put_in_T1(i)  = E[-(x_T1(i) + DeltaX_T1_T2 + m2)^+]
%
% The expectations are computed empirically using the simulated
% increments dx_T1_T2, but without looping over all inner samples
% for each outer scenario.

a_T1 = x_T1 + m2;

[call_in_T1, put_in_T1] = inner_call_put_from_sorted_increments( ...
    a_T1, dx_T1_T2);

%% Price selected exotic and MC standard error
%
% We explicitly build the Monte Carlo variable Y_i.
% Then:
%   price = discountExotic * mean(Y_i)
%   SE    = discountExotic * std(Y_i) / sqrt(N)

switch exoticType

    case 'CoC'
        % Y_i^CoC = [B(T1,T2) * C_T1^{(i)}(K2) - K1]^+
        Y_exotic = max(B12 .* call_in_T1 - K1, 0);
        discountExotic = B_T1;

    case 'PoP'
        % Y_i^PoP = [K1 - B(T1,T2) * P_T1^{(i)}(K2)]^+
        Y_exotic = max(K1 - B12 .* put_in_T1, 0);
        discountExotic = B_T1;

    case 'Chooser'
        % Y_i^Chooser = max(C_T1^{(i)}(K2), P_T1^{(i)}(K2))
        Y_exotic = max(call_in_T1, put_in_T1);
        discountExotic = B_T2;
end

N_outer = numel(Y_exotic);

price_exotic = discountExotic * mean(Y_exotic);

mcStdError = discountExotic * std(Y_exotic, 0) / sqrt(N_outer);

ci95HalfWidth = 1.96 * mcStdError;
ci95Low  = price_exotic - ci95HalfWidth;
ci95High = price_exotic + ci95HalfWidth;
%% Diagnostics

diagnostics = struct();

diagnostics.exoticType = exoticType;
diagnostics.K1 = K1;
diagnostics.K2 = K2;
diagnostics.F_T0_T2 = F_T0_T2;
diagnostics.m2 = m2;

diagnostics.N_outer = numel(x_T1);
diagnostics.N_inner = numel(dx_T1_T2);
diagnostics.N_sim = numel(x_T1);

diagnostics.method = 'sorted_inner_increments_no_grid';
diagnostics.num_grid_points = 0;
diagnostics.num_grid_points_input_ignored = num_grid_points;

diagnostics.Y_exotic_mean = mean(Y_exotic);
diagnostics.Y_exotic_std  = std(Y_exotic, 0);

diagnostics.discountExotic = discountExotic;

diagnostics.mcStdError = mcStdError;
diagnostics.ci95_half_width = ci95HalfWidth;
diagnostics.ci95_low = ci95Low;
diagnostics.ci95_high = ci95High;

diagnostics.priceMean = price_exotic;

diagnostics.call_in_T1_min = min(call_in_T1);
diagnostics.call_in_T1_max = max(call_in_T1);
diagnostics.put_in_T1_min = min(put_in_T1);
diagnostics.put_in_T1_max = max(put_in_T1);

diagnostics.mean_dx_T1_T2 = mean(dx_T1_T2);

% Empirical put-call parity check:
%
% Since call and put are computed with the same empirical distribution
% of dx_T1_T2, the following identity should hold numerically:
%
%   call_in_T1 - put_in_T1 = x_T1 + m2 + mean(dx_T1_T2)
%
% No zero-mean assumption is used here.
diagnostics.inner_parity_max_abs_error = max(abs( ...
    call_in_T1 - put_in_T1 - ...
    (x_T1 + m2 + diagnostics.mean_dx_T1_T2) ...
    ));

% Vanilla prices obtained by recycling the inner increment distribution.
% These are useful benchmark quantities for the inner layer.
diagnostics.vanilla_call_T2_recycled_MC = B_T2 * mean(call_in_T1);
diagnostics.vanilla_put_T2_recycled_MC  = B_T2 * mean(put_in_T1);

% Direct paired MC diagnostics are available only when the two vectors
% have the same length.
sameLengthSamples = numel(x_T1) == numel(dx_T1_T2);
diagnostics.same_length_samples = sameLengthSamples;

if sameLengthSamples
    X_T2_from_steps = x_T1 + dx_T1_T2;

    diagnostics.vanilla_call_T2_direct_MC = ...
        B_T2 * mean(max( X_T2_from_steps + m2, 0));

    diagnostics.vanilla_put_T2_direct_MC = ...
        B_T2 * mean(max(-X_T2_from_steps - m2, 0));

    % Same diagnostic as in the mesh version.
    % When m2 = 0 and mean(dx_T1_T2) is close to zero, this is the
    % standard chooser identity used as a sanity check.
    diagnostics.chooser_identity_MC = B_T2 * ( ...
        mean(max(X_T2_from_steps + m2, 0)) + ...
        mean(max(-x_T1 - m2, 0)) );
else
    diagnostics.vanilla_call_T2_direct_MC = NaN;
    diagnostics.vanilla_put_T2_direct_MC = NaN;
    diagnostics.chooser_identity_MC = NaN;
end

% Empirical chooser identity based on the computed inner values.
% This does not require mean(dx_T1_T2) = 0.
diagnostics.chooser_identity_empirical = B_T2 * mean( ...
    call_in_T1 + max(-(x_T1 + m2 + diagnostics.mean_dx_T1_T2), 0) ...
    );

diagnostics.chooser_identity_empirical_error = abs( ...
    B_T2 * mean(max(call_in_T1, put_in_T1)) - ...
    diagnostics.chooser_identity_empirical ...
    );

end


function [callMean, putMean] = inner_call_put_from_sorted_increments(a, z)
%INNER_CALL_PUT_FROM_SORTED_INCREMENTS
% Computes the empirical inner expectations:
%
%   callMean(i) = mean(max( a(i) + z, 0))
%   putMean(i)  = mean(max(-a(i) - z, 0))
%
% Inputs:
%   a : N_outer x 1 vector
%       In the project, a = x_T1 + m2.
%
%   z : N_inner x 1 vector
%       In the project, z = DeltaX_T1_T2.
%
% No zero-mean assumption is used.

a = a(:);
z = z(:);

M = numel(z);

% Unique sorted inner increments.
% The grouping is needed because discretize requires strictly
% increasing bin edges.
[zUnique, ~, groupIdx] = unique(z);

counts = accumarray(groupIdx, 1);
sums   = accumarray(groupIdx, z);

cumCount = [0; cumsum(counts)];
cumSum   = [0; cumsum(sums)];

totalSum = sum(z);

% Exercise threshold:
%
% call active when a + z > 0  <=>  z > -a
% put  active when a + z < 0  <=>  z < -a
threshold = -a;

edges = [-Inf; zUnique; Inf];

% idx(i) is the number of unique z-levels <= threshold(i).
idx = discretize(threshold, edges) - 1;

idx(isnan(idx)) = 0;
idx = max(0, min(idx, numel(zUnique)));

leftCount = cumCount(idx + 1);
leftSum   = cumSum(idx + 1);

rightCount = M - leftCount;
rightSum   = totalSum - leftSum;

% callMean(i) = mean((a(i) + z)^+) over z > -a(i)
callMean = (rightCount .* a + rightSum) ./ M;

% putMean(i) = mean((-a(i) - z)^+) over z <= -a(i)
putMean = (-leftCount .* a - leftSum) ./ M;

% Numerical protection only.
callMean = max(callMean, 0);
putMean  = max(putMean, 0);

end