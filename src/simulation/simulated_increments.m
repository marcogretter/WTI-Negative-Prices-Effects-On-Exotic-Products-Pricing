function [xtarget]= simulated_increments(u, invCDF)
% Returns the increment vector via interpolation in the set on which the
% invCDF is defined (MID) and via exponential approximation in the tails.

xtarget = zeros(size(u));

idxLeft  = u <= invCDF.Pb;
idxRight = u >= invCDF.Pe;
idxMid   = ~idxLeft & ~idxRight;

% Left tail
xtarget(idxLeft) = invCDF.xBlock(1) + (1 / invCDF.lambdaMinus) .* log(u(idxLeft) ./ invCDF.Pb);

% Right tail
xtarget(idxRight) = invCDF.xBlock(end) - (1 / invCDF.lambdaPlus) .* log((1 - u(idxRight)) ./ (1 - invCDF.Pe));

% Middle spline interpolation
xtarget(idxMid) = ppval(spline(invCDF.PBlock, invCDF.xBlock), u(idxMid));