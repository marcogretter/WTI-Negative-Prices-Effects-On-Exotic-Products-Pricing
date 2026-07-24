function point4 = run_point4_exotic_pricing(simIncrements, vol, point4Params)
% Price CoC, PoP and Chooser for AB/GL/MA using grid and no-grid methods.
% Benchmarks timing, computes speedup, and captures inner parity errors.
% Also runs the MA analytical check via check_price_exotics_MA.
%
% Inputs:
%   simIncrements  struct from build_simIncrements
%   vol            struct from point2_calibrate_vol_surface
%   point4Params   struct: num_grid_points (def 2000), K1_CoC (def 0.3), K1_PoP (def 0.3)
%
% Output struct point4:
%   tables.gridPrices       model x exotic prices (stochastic mesh)
%   tables.noGridPrices     model x exotic prices (no-grid)
%   tables.timingComparison timing and speedup aggregated per model
%   diagnostics             per-model/exotic structs, parity errors, MA analytical check

if isfield(point4Params, 'num_grid_points')
    num_grid_points = point4Params.num_grid_points;
else
    num_grid_points = 2000;
end

if isfield(point4Params, 'K1_CoC')
    K1_CoC = point4Params.K1_CoC;
else
    K1_CoC = 0.3;
end

if isfield(point4Params, 'K1_PoP')
    K1_PoP = point4Params.K1_PoP;
else
    K1_PoP = 0.3;
end

B_T1    = simIncrements.common.B_T1;
B_T2    = simIncrements.common.B_T2;
F_T0_T2 = simIncrements.common.F_T0_T2;
tGrid_d = simIncrements.common.tGrid_d;
K2      = F_T0_T2;

modelNames  = {'AB', 'GL', 'MA'};
exoticNames = {'CoC', 'PoP', 'Chooser'};
nModels     = numel(modelNames);
nExotics    = numel(exoticNames);

gridPricesMat   = zeros(nModels, nExotics);
noGridPricesMat = zeros(nModels, nExotics);
gridTimeMat     = zeros(nModels, nExotics);
noGridTimeMat   = zeros(nModels, nExotics);
parityErrorMat  = nan(nModels,   nExotics);
% Absolute MC errors
gridSEMat   = nan(nModels, nExotics);
noGridSEMat = nan(nModels, nExotics);

gridCI95HalfWidthMat = nan(nModels, nExotics);
noGridCI95HalfWidthMat = nan(nModels, nExotics);

% Forward-based bps
gridSEBpsForwardMat = nan(nModels, nExotics);
noGridSEBpsForwardMat = nan(nModels, nExotics);

gridCI95HalfWidthBpsForwardMat = nan(nModels, nExotics);
noGridCI95HalfWidthBpsForwardMat = nan(nModels, nExotics);

gridPass5bpsForwardMat = false(nModels, nExotics);
noGridPass5bpsForwardMat = false(nModels, nExotics);

gridPass95_5bpsForwardMat = false(nModels, nExotics);
noGridPass95_5bpsForwardMat = false(nModels, nExotics);

% Price-based bps
gridSEBpsPriceMat = nan(nModels, nExotics);
noGridSEBpsPriceMat = nan(nModels, nExotics);

gridCI95HalfWidthBpsPriceMat = nan(nModels, nExotics);
noGridCI95HalfWidthBpsPriceMat = nan(nModels, nExotics);

gridPass5bpsPriceMat = false(nModels, nExotics);
noGridPass5bpsPriceMat = false(nModels, nExotics);

gridPass95_5bpsPriceMat = false(nModels, nExotics);
noGridPass95_5bpsPriceMat = false(nModels, nExotics);

% Forward-based bps reference
forwardBpsBase = abs(F_T0_T2);

if forwardBpsBase <= 0 || ~isfinite(forwardBpsBase)
    error('Invalid bps base: F_T0_T2 must be positive and finite.');
end

% Price-based bps protection
priceFloor = 1e-12;

% Threshold and 95% Gaussian confidence interval multiplier
thresholdBps = 5;
z95 = 1.96;
point4.diagnostics = struct();

for i = 1:nModels
    modelName = modelNames{i};
    x_T0_T1   = simIncrements.(modelName).increments.T0_T1;
    x_T1_T2   = simIncrements.(modelName).increments.T1_T2;

    for j = 1:nExotics
        exoticType = exoticNames{j};
        switch exoticType
            case 'CoC',     K1 = K1_CoC;
            case 'PoP',     K1 = K1_PoP;
            case 'Chooser', K1 = [];
        end

        t0 = tic;
        [priceGrid, diagGrid] = price_exotic_mesh_from_increments( ...
            exoticType, x_T0_T1, x_T1_T2, ...
            B_T1, B_T2, F_T0_T2, K2, K1, num_grid_points);
        gridTimeMat(i, j)   = toc(t0);
        gridPricesMat(i, j) = priceGrid;
        gridSEMat(i, j) = diagGrid.mcStdError;

        gridSEMat(i, j) = diagGrid.mcStdError;
        gridCI95HalfWidthMat(i, j) = z95 * diagGrid.mcStdError;
        
        priceBaseGrid = max(abs(priceGrid), priceFloor);
        
        % Forward-based bps
        gridSEBpsForwardMat(i, j) = ...
            1e4 * diagGrid.mcStdError / forwardBpsBase;
        
        gridCI95HalfWidthBpsForwardMat(i, j) = ...
            1e4 * gridCI95HalfWidthMat(i, j) / forwardBpsBase;
        
        gridPass5bpsForwardMat(i, j) = ...
            gridSEBpsForwardMat(i, j) < thresholdBps;
        
        gridPass95_5bpsForwardMat(i, j) = ...
            gridCI95HalfWidthBpsForwardMat(i, j) < thresholdBps;
        
        % Price-based bps
        gridSEBpsPriceMat(i, j) = ...
            1e4 * diagGrid.mcStdError / priceBaseGrid;
        
        gridCI95HalfWidthBpsPriceMat(i, j) = ...
            1e4 * gridCI95HalfWidthMat(i, j) / priceBaseGrid;
        
        gridPass5bpsPriceMat(i, j) = ...
            gridSEBpsPriceMat(i, j) < thresholdBps;
        
        gridPass95_5bpsPriceMat(i, j) = ...
            gridCI95HalfWidthBpsPriceMat(i, j) < thresholdBps;

        t0 = tic;
        [priceNoGrid, diagNoGrid] = price_exotic_no_grid_from_increments( ...
            exoticType, x_T0_T1, x_T1_T2, ...
            B_T1, B_T2, F_T0_T2, K2, K1);
        noGridTimeMat(i, j)   = toc(t0);
        noGridPricesMat(i, j) = priceNoGrid;
        noGridSEMat(i, j) = diagNoGrid.mcStdError;
        noGridCI95HalfWidthMat(i, j) = z95 * diagNoGrid.mcStdError;
        
        priceBaseNoGrid = max(abs(priceNoGrid), priceFloor);
        
        % Forward-based bps
        noGridSEBpsForwardMat(i, j) = ...
            1e4 * diagNoGrid.mcStdError / forwardBpsBase;
        
        noGridCI95HalfWidthBpsForwardMat(i, j) = ...
            1e4 * noGridCI95HalfWidthMat(i, j) / forwardBpsBase;
        
        noGridPass5bpsForwardMat(i, j) = ...
            noGridSEBpsForwardMat(i, j) < thresholdBps;
        
        noGridPass95_5bpsForwardMat(i, j) = ...
            noGridCI95HalfWidthBpsForwardMat(i, j) < thresholdBps;
        
        % Price-based bps
        noGridSEBpsPriceMat(i, j) = ...
            1e4 * diagNoGrid.mcStdError / priceBaseNoGrid;
        
        noGridCI95HalfWidthBpsPriceMat(i, j) = ...
            1e4 * noGridCI95HalfWidthMat(i, j) / priceBaseNoGrid;
        
        noGridPass5bpsPriceMat(i, j) = ...
            noGridSEBpsPriceMat(i, j) < thresholdBps;
        
        noGridPass95_5bpsPriceMat(i, j) = ...
            noGridCI95HalfWidthBpsPriceMat(i, j) < thresholdBps;
        
        point4.diagnostics.(modelName).(exoticType).grid   = diagGrid;
        point4.diagnostics.(modelName).(exoticType).noGrid = diagNoGrid;

        if isfield(diagNoGrid, 'inner_parity_max_abs_error')
            parityErrorMat(i, j) = diagNoGrid.inner_parity_max_abs_error;
        end
    end
end

% MA analytical check
alpha_MA    = vol.MA.alpha;
beta_MA     = vol.MA.beta;
sigma_T1_MA = simIncrements.MA.sigma.T1;
sigma_T2_MA = simIncrements.MA.sigma.T2;

for j = 1:nExotics
    exoticType = exoticNames{j};
    [priceNum, priceAna] = check_price_exotics_MA( ...
        alpha_MA, beta_MA, tGrid_d, sigma_T1_MA, sigma_T2_MA, ...
        B_T1, B_T2, ...
        simIncrements.MA.increments.T0_T1, simIncrements.MA.increments.T1_T2, ...
        exoticType);
    point4.diagnostics.MA_analytical.(exoticType).numerical = priceNum;
    point4.diagnostics.MA_analytical.(exoticType).analytic  = priceAna;
    point4.diagnostics.MA_analytical.(exoticType).absDiff   = abs(priceNum - priceAna);
end

% Price tables

point4.tables.gridPrices = table( ...
    modelNames(:), ...
    gridPricesMat(:, 1), gridPricesMat(:, 2), gridPricesMat(:, 3), ...
    'VariableNames', {'Model', 'CoC', 'PoP', 'Chooser'});

point4.tables.noGridPrices = table( ...
    modelNames(:), ...
    noGridPricesMat(:, 1), noGridPricesMat(:, 2), noGridPricesMat(:, 3), ...
    'VariableNames', {'Model', 'CoC', 'PoP', 'Chooser'});
% Grid vs no-grid price differences

priceDiffMat = noGridPricesMat - gridPricesMat;
absPriceDiffMat = abs(priceDiffMat);

% Difference in bps with respect to the 12M forward
priceDiffBpsForwardMat = 1e4 * priceDiffMat ./ forwardBpsBase;
absPriceDiffBpsForwardMat = abs(priceDiffBpsForwardMat);

% Difference in bps with respect to the grid price
priceDiffBpsPriceMat = 1e4 * priceDiffMat ./ max(abs(gridPricesMat), priceFloor);
absPriceDiffBpsPriceMat = abs(priceDiffBpsPriceMat);

point4.tables.gridNoGridPriceDiff = table( ...
    modelNames(:), ...
    priceDiffMat(:, 1), priceDiffMat(:, 2), priceDiffMat(:, 3), ...
    absPriceDiffMat(:, 1), absPriceDiffMat(:, 2), absPriceDiffMat(:, 3), ...
    absPriceDiffBpsForwardMat(:, 1), absPriceDiffBpsForwardMat(:, 2), absPriceDiffBpsForwardMat(:, 3), ...
    absPriceDiffBpsPriceMat(:, 1), absPriceDiffBpsPriceMat(:, 2), absPriceDiffBpsPriceMat(:, 3), ...
    'VariableNames', { ...
        'Model', ...
        'CoC_NoGridMinusGrid', 'PoP_NoGridMinusGrid', 'Chooser_NoGridMinusGrid', ...
        'CoC_AbsDiff', 'PoP_AbsDiff', 'Chooser_AbsDiff', ...
        'CoC_AbsDiff_bps_forward', 'PoP_AbsDiff_bps_forward', 'Chooser_AbsDiff_bps_forward', ...
        'CoC_AbsDiff_bps_price', 'PoP_AbsDiff_bps_price', 'Chooser_AbsDiff_bps_price' ...
    });

point4.tables.gridMCErrorForward = table( ...
    modelNames(:), ...
    gridSEMat(:, 1), gridSEMat(:, 2), gridSEMat(:, 3), ...
    gridSEBpsForwardMat(:, 1), gridSEBpsForwardMat(:, 2), gridSEBpsForwardMat(:, 3), ...
    gridCI95HalfWidthBpsForwardMat(:, 1), gridCI95HalfWidthBpsForwardMat(:, 2), gridCI95HalfWidthBpsForwardMat(:, 3), ...
    gridPass5bpsForwardMat(:, 1), gridPass5bpsForwardMat(:, 2), gridPass5bpsForwardMat(:, 3), ...
    gridPass95_5bpsForwardMat(:, 1), gridPass95_5bpsForwardMat(:, 2), gridPass95_5bpsForwardMat(:, 3), ...
    'VariableNames', { ...
        'Model', ...
        'CoC_SE', 'PoP_SE', 'Chooser_SE', ...
        'CoC_SE_bps_forward', 'PoP_SE_bps_forward', 'Chooser_SE_bps_forward', ...
        'CoC_CI95HalfWidth_bps_forward', 'PoP_CI95HalfWidth_bps_forward', 'Chooser_CI95HalfWidth_bps_forward', ...
        'CoC_Pass_5bps_forward', 'PoP_Pass_5bps_forward', 'Chooser_Pass_5bps_forward', ...
        'CoC_Pass95_5bps_forward', 'PoP_Pass95_5bps_forward', 'Chooser_Pass95_5bps_forward' ...
    });

point4.tables.noGridMCErrorForward = table( ...
    modelNames(:), ...
    noGridSEMat(:, 1), noGridSEMat(:, 2), noGridSEMat(:, 3), ...
    noGridSEBpsForwardMat(:, 1), noGridSEBpsForwardMat(:, 2), noGridSEBpsForwardMat(:, 3), ...
    noGridCI95HalfWidthBpsForwardMat(:, 1), noGridCI95HalfWidthBpsForwardMat(:, 2), noGridCI95HalfWidthBpsForwardMat(:, 3), ...
    noGridPass5bpsForwardMat(:, 1), noGridPass5bpsForwardMat(:, 2), noGridPass5bpsForwardMat(:, 3), ...
    noGridPass95_5bpsForwardMat(:, 1), noGridPass95_5bpsForwardMat(:, 2), noGridPass95_5bpsForwardMat(:, 3), ...
    'VariableNames', { ...
        'Model', ...
        'CoC_SE', 'PoP_SE', 'Chooser_SE', ...
        'CoC_SE_bps_forward', 'PoP_SE_bps_forward', 'Chooser_SE_bps_forward', ...
        'CoC_CI95HalfWidth_bps_forward', 'PoP_CI95HalfWidth_bps_forward', 'Chooser_CI95HalfWidth_bps_forward', ...
        'CoC_Pass_5bps_forward', 'PoP_Pass_5bps_forward', 'Chooser_Pass_5bps_forward', ...
        'CoC_Pass95_5bps_forward', 'PoP_Pass95_5bps_forward', 'Chooser_Pass95_5bps_forward' ...
    });
point4.tables.gridMCErrorPrice = table( ...
    modelNames(:), ...
    gridSEMat(:, 1), gridSEMat(:, 2), gridSEMat(:, 3), ...
    gridSEBpsPriceMat(:, 1), gridSEBpsPriceMat(:, 2), gridSEBpsPriceMat(:, 3), ...
    gridCI95HalfWidthBpsPriceMat(:, 1), gridCI95HalfWidthBpsPriceMat(:, 2), gridCI95HalfWidthBpsPriceMat(:, 3), ...
    gridPass5bpsPriceMat(:, 1), gridPass5bpsPriceMat(:, 2), gridPass5bpsPriceMat(:, 3), ...
    gridPass95_5bpsPriceMat(:, 1), gridPass95_5bpsPriceMat(:, 2), gridPass95_5bpsPriceMat(:, 3), ...
    'VariableNames', { ...
        'Model', ...
        'CoC_SE', 'PoP_SE', 'Chooser_SE', ...
        'CoC_SE_bps_price', 'PoP_SE_bps_price', 'Chooser_SE_bps_price', ...
        'CoC_CI95HalfWidth_bps_price', 'PoP_CI95HalfWidth_bps_price', 'Chooser_CI95HalfWidth_bps_price', ...
        'CoC_Pass_5bps_price', 'PoP_Pass_5bps_price', 'Chooser_Pass_5bps_price', ...
        'CoC_Pass95_5bps_price', 'PoP_Pass95_5bps_price', 'Chooser_Pass95_5bps_price' ...
    });

point4.tables.noGridMCErrorPrice = table( ...
    modelNames(:), ...
    noGridSEMat(:, 1), noGridSEMat(:, 2), noGridSEMat(:, 3), ...
    noGridSEBpsPriceMat(:, 1), noGridSEBpsPriceMat(:, 2), noGridSEBpsPriceMat(:, 3), ...
    noGridCI95HalfWidthBpsPriceMat(:, 1), noGridCI95HalfWidthBpsPriceMat(:, 2), noGridCI95HalfWidthBpsPriceMat(:, 3), ...
    noGridPass5bpsPriceMat(:, 1), noGridPass5bpsPriceMat(:, 2), noGridPass5bpsPriceMat(:, 3), ...
    noGridPass95_5bpsPriceMat(:, 1), noGridPass95_5bpsPriceMat(:, 2), noGridPass95_5bpsPriceMat(:, 3), ...
    'VariableNames', { ...
        'Model', ...
        'CoC_SE', 'PoP_SE', 'Chooser_SE', ...
        'CoC_SE_bps_price', 'PoP_SE_bps_price', 'Chooser_SE_bps_price', ...
        'CoC_CI95HalfWidth_bps_price', 'PoP_CI95HalfWidth_bps_price', 'Chooser_CI95HalfWidth_bps_price', ...
        'CoC_Pass_5bps_price', 'PoP_Pass_5bps_price', 'Chooser_Pass_5bps_price', ...
        'CoC_Pass95_5bps_price', 'PoP_Pass95_5bps_price', 'Chooser_Pass95_5bps_price' ...
    });
% Timing comparison aggregated over exotics per model

totalGridTime   = sum(gridTimeMat,   2);
totalNoGridTime = sum(noGridTimeMat, 2);
speedUp         = totalGridTime ./ totalNoGridTime;

point4.tables.timingComparison = table( ...
    modelNames(:), ...
    totalGridTime, totalNoGridTime, speedUp, ...
    'VariableNames', {'Model', 'GridTime_s', 'NoGridTime_s', 'SpeedUp'});

% Inner parity error summary

point4.diagnostics.innerParityErrors = table( ...
    modelNames(:), ...
    parityErrorMat(:, 1), parityErrorMat(:, 2), parityErrorMat(:, 3), ...
    'VariableNames', {'Model', 'CoC', 'PoP', 'Chooser'});

% Backward-compatible aliases

point4.tables.gridMCError   = point4.tables.gridMCErrorForward;
point4.tables.noGridMCError = point4.tables.noGridMCErrorForward;

disp('Point 4: exotic prices — grid method');
disp(point4.tables.gridPrices);

disp('Point 4: exotic prices — no-grid method');
disp(point4.tables.noGridPrices);

disp('Point 4: timing comparison');
disp(point4.tables.timingComparison);

disp('Point 4: MC standard errors — grid method, forward-based bps');
disp(point4.tables.gridMCErrorForward);

disp('Point 4: MC standard errors — no-grid method, forward-based bps');
disp(point4.tables.noGridMCErrorForward);

disp('Point 4: MC standard errors — grid method, price-based bps');
disp(point4.tables.gridMCErrorPrice);

disp('Point 4: MC standard errors — no-grid method, price-based bps');
disp(point4.tables.noGridMCErrorPrice);

disp('Point 4: grid vs no-grid price differences');
disp(point4.tables.gridNoGridPriceDiff);
end
