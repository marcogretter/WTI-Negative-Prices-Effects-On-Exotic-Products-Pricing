function [zero_rates]=plot_calibration_df_and_forwards(B_curve, F_curve, futureExpiries, valueDate)
    if nargin < 4
        valueDate = datetime(2020, 6, 2);
    end

    B_curve = B_curve(:);
    F_curve = F_curve(:);
    futureExpiries = futureExpiries(:);

    TTM_all = yearfrac(valueDate, futureExpiries, 3);

    % evaluation of the filter
    valid_idx = isfinite(B_curve) & B_curve > 0 & isfinite(F_curve) & TTM_all > 0;
    
    % apply the filter to get B and F
    B_plot = B_curve(valid_idx);
    F_plot = F_curve(valid_idx);
    
    % Evaluation of the TTM with yearfrac ACT/365
    TTM = TTM_all(valid_idx);
    
    
    % Evaluation of the Zero Rates (continuous compounding)
    zero_rates = -log(B_plot) ./ TTM;
    
    % Plot
    figure('Name', 'Calibration Discount, Forward and Zero Rates', 'NumberTitle', 'off');
    
    % Plot Discount Factors
    subplot(3,1,1);
    plot(TTM, B_plot, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
    grid on; 
    xlabel('Time to Maturity (years)');
    ylabel('Discount Factor');
    title('Calibration Discount Factors (B)');
    
    % Plot Forward Prices
    subplot(3,1,2);
    plot(TTM, F_plot, '-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0.85 0.33 0.1]);
    grid on; 
    xlabel('Time to Maturity (years)');
    ylabel('Forward Price');
    title('Calibration Forward Prices (F)');
    
    % Plot Zero Rates
    subplot(3,1,3);
    % Multiply by 100 to get rates in %
    plot(TTM, zero_rates * 100, '-^', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0.47 0.67 0.19]);
    grid on; 
    xlabel('Time to Maturity (years)');
    ylabel('Zero Rate (%)');
    title('Calibration Zero Rates');
end
