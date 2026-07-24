function invCDF = build_inverse_cdf_spline(x_grid, P_grid, lambdaMinusTheory, lambdaPlusTheory)
%BUILD_INVERSE_CDF_SPLINE Build an inverse CDF interpolant from a numerical CDF.
%
% The function selects the largest valid monotone block of the reconstructed
% CDF, builds a spline inverse on that block, and stores exponential tail
% coefficients for inverse CDF extrapolation.
%
% If theoretical tail coefficients are provided, they are used directly.
% Otherwise, the coefficients are estimated from the boundary points of the
% valid CDF block.

x = x_grid(:);
P = P_grid(:);

[x, ord] = sort(x);
P = P(ord);

epsCDF = 1e-12;
tolInc = 1e-10;

validPoint = isfinite(x) & isfinite(P) & P > epsCDF & P < 1 - epsCDF;

[istart, iend] = largestIncreasingBlock(P, validPoint, tolInc);

if isempty(istart) || isempty(iend)
    error('No valid block found for inverse CDF construction.');
end

xBlock = x(istart:iend);
PBlock = P(istart:iend);

keep = [true; diff(PBlock) > tolInc];

xBlock = xBlock(keep);
PBlock = PBlock(keep);

if length(xBlock) < 4
    error('Not enough points in the valid block for spline construction.');
end

xb = xBlock(1);
xe = xBlock(end);

Pb = PBlock(1);
Pe = PBlock(end);

% Exponential extrapolation parameters for left and right tails
%
% If theoretical tail coefficients are provided, use them.
% Otherwise estimate them from the first/last valid CDF points.

if nargin >= 4 && ~isempty(lambdaMinusTheory) && ~isempty(lambdaPlusTheory)

    lambdaMinus = lambdaMinusTheory;
    lambdaPlus  = lambdaPlusTheory;

    tailMode = 'theoretical';

else

    lambdaMinus = (log(PBlock(2)) - log(PBlock(1))) / ...
        (xBlock(2) - xBlock(1));

    lambdaPlus = (log(PBlock(end)) - log(PBlock(end-1))) / ...
        (xBlock(end) - xBlock(end-1));

    tailMode = 'estimated';

end

if lambdaMinus <= 0 || lambdaPlus <= 0
    warning('Non-positive tail coefficient detected: lambdaMinus = %.6g, lambdaPlus = %.6g.', ...
        lambdaMinus, lambdaPlus);
end

ppInv = spline(PBlock, xBlock);

invCDF.ppInv = ppInv;

invCDF.xb = xb;
invCDF.xe = xe;

invCDF.Pb = Pb;
invCDF.Pe = Pe;

invCDF.lambdaMinus = lambdaMinus;
invCDF.lambdaPlus  = lambdaPlus;

invCDF.tailMode = tailMode;

invCDF.xBlock = xBlock;
invCDF.PBlock = PBlock;
end


function [bestStart, bestEnd] = largestIncreasingBlock(P, validPoint, tolInc)
% finds the largest consecutive block [xb, xe] such that:
% validPoint(k) == true for every point in the block
% and P is increasing by at least tolInc across the block

P = P(:);
validPoint = validPoint(:);

N = length(P);

if N == 0
    bestStart = [];
    bestEnd = [];
    return;
end

goodedge = validPoint(1:end-1) & validPoint(2:end) & (diff(P) > tolInc);

blockstartMask = validPoint & [true; ~goodedge];

blockendMask = validPoint & [~goodedge; true];

starts = find(blockstartMask);
ends   = find(blockendMask);

if isempty(starts) || isempty(ends)
    bestStart = [];
    bestEnd = [];
    return;
end

lengths = ends - starts + 1;

[~, idxBest] = max(lengths);

bestStart = starts(idxBest);
bestEnd = ends(idxBest);
end