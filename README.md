# WTI Exotic Option Pricing with Linear Additive Models in MATLAB

This repository contains a MATLAB library for the calibration, simulation,
pricing, and risk management of WTI crude oil derivatives under three linear
additive models:

- Additive Bachelier;
- Minimal Additive;
- Generalized Logistic.

The project implements a complete quantitative derivatives workflow:

1. extraction of synthetic forwards and discount factors from option prices;
2. calibration of the WTI implied-volatility surface;
3. simulation of additive processes through Lewis-FFT inversion;
4. pricing of forward-start and exotic options;
5. derivation of analytical benchmarks for the Minimal Additive model;
6. construction and backtesting of a Delta–Vega hedge.

The market data refer to WTI options observed on 2 June 2020, during the period
of extreme oil-market volatility following the Covid-19 shock.

## Project Motivation

The oil market in 2020 required pricing models capable of handling very low and
potentially negative forward prices.

Traditional lognormal models are not naturally suited to this environment.

The project therefore uses a linear additive representation:

```text
forward_price(t, T) =
    initial_forward_price(0, T)
    + additive_factor_t
```

Unlike a multiplicative lognormal model, this specification allows the forward
price to take any real value.

The additive structure also provides:

- independent increments;
- explicit characteristic functions;
- flexible asymmetric distributions;
- exponential tails;
- Fourier-based pricing and simulation;
- tractable analytical results in selected cases.

## Main Objectives

The project investigates the following questions:

- How can discount factors and forward prices be extracted directly from WTI
  option quotes?
- Can different linear additive models reproduce the observed volatility smile?
- How can additive increments be simulated efficiently from their characteristic
  functions?
- How accurate is Lewis-FFT Monte Carlo pricing?
- How sensitive are exotic prices to the selected distributional model?
- Can analytical formulas be derived for the Minimal Additive model?
- Can the exotic portfolio be effectively hedged with vanilla options and
  futures?

## Data

The market dataset contains:

```text
Data/
|
|-- datacalls/
|   `-- one file for each option expiry
|
|-- dataputs/
|   `-- one file for each option expiry
|
`-- Expiries_Futures
```

The call and put files contain option mid prices across a grid of strikes and
valuation dates.

The futures-expiry file maps each option maturity to the corresponding WTI
futures contract.

The principal valuation date is:

```text
2 June 2020
```

Additional market dates are used in the hedging backtest:

```text
9 June 2020

16 June 2020
```

## Models

The project compares three linear additive models.

### Additive Bachelier Model

The Additive Bachelier model introduces asymmetry and non-Gaussian tails through
a normal mean-variance mixture.

Its main shape parameters are:

```text
eta =
    skew parameter

k =
    volatility-of-volatility parameter
```

The benchmark subordinator parameter is fixed at:

```text
alpha = 0.5
```

The model is able to produce:

- asymmetric distributions;
- excess kurtosis;
- smooth implied-volatility smiles;
- real-valued forward prices.

### Minimal Additive Model

The normalized Minimal Additive distribution is asymmetric Laplace.

Its shape is controlled by:

```text
alpha_MA =
    left-tail coefficient

beta_MA =
    right-tail coefficient
```

The distribution has:

- exponential left and right tails;
- a sharp central cusp;
- elementary European option formulas;
- a probability atom in future increments;
- analytically invertible increment distributions.

The last property makes the model especially useful as a benchmark for the
FFT-based simulation engine.

### Generalized Logistic Model

The Generalized Logistic model uses a Type-IV generalized logistic
distribution.

Its shape parameters are:

```text
alpha_GL

beta_GL
```

The model produces:

- asymmetric exponential tails;
- a smooth central density;
- flexible intermediate-tail behavior;
- a characteristic function expressed through complex Gamma functions.

A dedicated complex log-Gamma implementation based on the Lanczos approximation
is included because the standard MATLAB `gammaln` function does not support
complex arguments.

## Linear Additive Representation

For all three models, the stochastic component is represented as:

```text
additive_factor_t =
    sigma_t
    * sqrt(t)
    * normalized_factor_t
```

The parameter:

```text
sigma_t
```

controls the maturity-dependent scale.

The distribution of:

```text
normalized_factor_t
```

determines the shape of the smile and the asymmetry of the tails.

This separation between scale and shape is the basis of the cascade calibration
procedure.

## Part 0 — Asymptotic Tail Analysis

The first section derives the exponential decay coefficients of the three model
densities.

For a generic additive factor, the tails behave approximately as:

```text
density(x)
    proportional to
exp(left_tail_coefficient * x)

for x tending to minus infinity
```

and:

```text
density(x)
    proportional to
exp(-right_tail_coefficient * x)

for x tending to plus infinity
```

The tail coefficients also define the analyticity strip of the characteristic
function.

For the scaled additive factor:

```text
time_dependent_tail_coefficient =
    normalized_tail_coefficient
    / (
        sigma_t
        * sqrt(t)
    )
```

The three model densities are compared under common normalized tail
coefficients:

```text
left_tail_coefficient = 0.9

right_tail_coefficient = 1.5
```

A log-linear density plot highlights the different behavior of the models near
the center and in the intermediate tails.

Main qualitative differences include:

- the sharp cusp of the Minimal Additive density;
- the smooth central region of the Generalized Logistic density;
- the smooth but more pronounced intermediate decay of the Additive Bachelier
  density.

## Part 1 — Synthetic Forward Calibration

Forward prices and discount factors are extracted from European call and put
prices using put-call parity.

For a fixed maturity and strike:

```text
call_price
    - put_price
    =
discount_factor
    * (
        forward_price
        - strike
    )
```

Defining:

```text
synthetic_forward_value =
    call_price
    - put_price
```

produces a linear relationship in the strike:

```text
synthetic_forward_value =
    regression_intercept
    + regression_slope
      * strike
```

The regression coefficients imply:

```text
discount_factor =
    -regression_slope
```

and:

```text
forward_price =
    regression_intercept
    / discount_factor
```

## Liquidity Filtering

Before estimating the synthetic-forward regression:

- calls and puts are matched on common strikes;
- options with mid prices below 0.1 are removed;
- maturities with fewer than three common liquid strikes are discarded;
- expired contracts are excluded.

The quality of each regression is evaluated through its coefficient of
determination.

## Calibrated Forward Curve

The market-implied forward prices obtained in the project are approximately:

```text
37.3505
37.9091
38.4689
39.0421
39.5549
40.0802
41.0900
42.1800
```

These values are used as model-independent inputs in all subsequent
calibrations.

## Discount Factors and Zero Rates

The calibrated discount factors are close to one because interest rates were
very low on the valuation date.

Continuously compounded zero rates are calculated as:

```text
zero_rate =
    -log(discount_factor)
    / time_to_maturity
```

The far end of the calibrated curve is interpreted carefully because fewer
liquid common strikes are available at longer maturities.

A high regression R-squared based on very few observations does not necessarily
imply a more reliable calibration.

## Absolute Dividends

Starting from consecutive forward prices and discount factors, the project
extracts deterministic absolute dividends.

A representative relation is:

```text
absolute_dividend =
    previous_discount_factor
    / current_discount_factor
    * previous_forward
    - current_forward
```

The corresponding intensity is:

```text
absolute_dividend_rate =
    absolute_dividend
    / maturity_interval
```

The resulting zero-rate and absolute-dividend curves are treated as deterministic
inputs in the rest of the project.

## Part 2 — Implied-Volatility Surface Calibration

The calibration is based on out-of-the-money options.

For each maturity:

- OTM calls are retained above the forward;
- OTM puts are retained below the forward;
- put prices are converted into call-equivalent prices;
- Bachelier implied volatilities are calculated;
- normalized moneyness is constructed.

## ATM Volatility Extraction

For each maturity, ATM normal volatility is estimated from the seven quotes
closest to zero moneyness.

A quadratic polynomial is fitted to:

```text
moneyness

implied_normal_volatility
```

The fitted polynomial is evaluated at zero moneyness.

When quotes are available on only one side of the forward, the nearest implied
volatility is used.

## Normalized Moneyness

Normalized moneyness is calculated as:

```text
normalized_moneyness =
    strike_minus_forward
    / (
        ATM_volatility
        * sqrt(time_to_maturity)
    )
```

This transformation allows options with different maturities and volatility
levels to be compared on a common scale.

## Cascade Calibration

The calibration follows three sequential stages:

```text
Stage 1:
    calibrate forwards and discount factors

Stage 2:
    extract ATM volatility term structure

Stage 3:
    calibrate distributional shape parameters
```

The shape parameters are calibrated globally across the complete OTM surface.

This is possible because the normalized pricing function is independent of
maturity once the scale has been separated from the distributional shape.

## Additive Bachelier Calibration

The Additive Bachelier parameters are calibrated by minimizing the squared
difference between market and model call-equivalent prices.

The initial guess is:

```text
eta_initial = 0.2

k_initial = 1.0
```

The calibrated values are:

```text
eta = -0.0624

k = 1.0180
```

The small negative value of `eta` indicates a mild left skew.

The value of `k` indicates a moderate degree of smile convexity.

## Minimal Additive Calibration

The Minimal Additive calibration identifies the relative left-right asymmetry
more reliably than the absolute scale of the two tail coefficients.

Multiple initial points produce nearly identical objective values and almost
identical ratios:

```text
alpha_MA / beta_MA
    approximately equals
1.026
```

The selected representative calibration is:

```text
alpha_MA = 1.0239

beta_MA = 0.9980
```

The dependence of the individual parameters on the optimizer starting point
reveals a numerical identifiability issue.

The available OTM quotes determine the relative tail asymmetry, but contain
limited information about a common rescaling of both coefficients.

## Generalized Logistic Calibration

The calibrated Generalized Logistic parameters are:

```text
alpha_GL = 1.0180

beta_GL = 0.3679
```

The pricing engine evaluates the characteristic function using the complex
log-Gamma routine implemented through the Lanczos approximation.

## Calibration Results

The global root mean squared pricing errors across 596 OTM options are:

```text
Additive Bachelier RMSE = 0.1007

Generalized Logistic RMSE = 0.1052

Minimal Additive RMSE = 0.1336
```

The Additive Bachelier and Generalized Logistic models provide similar global
fit quality.

The Minimal Additive model performs less accurately near the center of the
distribution because its asymmetric Laplace density has a cusp and no smooth
transition between the body and the exponential tails.

All three models perform less accurately at the longest maturity, where option
liquidity is lower.

## Part 3 — Increment Simulation

The additive processes have independent increments.

For two dates `s` and `t`, with `t > s`, the increment characteristic function is:

```text
increment_characteristic_function(s, t, u) =
    marginal_characteristic_function(t, u)
    / marginal_characteristic_function(s, u)
```

In logarithmic form:

```text
log_increment_CF =
    log_marginal_CF_at_t
    - log_marginal_CF_at_s
```

This property allows paths to be simulated sequentially between arbitrary reset
dates.

## Lewis-FFT CDF Reconstruction

Increment distributions are simulated by reconstructing their cumulative
distribution functions through Fourier inversion.

The procedure is:

1. evaluate the increment characteristic function on a Fourier grid;
2. shift the integration contour inside the analyticity strip;
3. apply the Fast Fourier Transform;
4. reconstruct the CDF on a finite spatial grid;
5. retain the monotone section of the numerical CDF;
6. interpolate the valid section;
7. extrapolate the tails;
8. invert the resulting CDF using uniform random variables.

## Complex Shift

The contour shift is selected inside the analyticity strip.

The implementation uses a value approximately halfway between zero and the
right boundary:

```text
complex_shift =
    right_tail_boundary
    / 2
```

This helps balance truncation and discretization errors.

## FFT Grid

The principal FFT configuration is:

```text
number_of_points = 32768

left_grid_boundary = -300

right_grid_boundary = 300
```

The range extends well beyond the option moneyness levels observed in the
dataset.

## Tail Extrapolation

The raw FFT CDF may not reach exactly zero and one at the grid boundaries.

An initial implementation estimated tail-decay rates from the boundary values of
the reconstructed CDF.

This method proved unstable, especially in the right tail where the CDF was
already numerically indistinguishable from one.

The final implementation uses the theoretical tail coefficients implied by the
characteristic function.

For the left tail:

```text
CDF(x) =
    CDF(left_boundary)
    * exp(
        left_tail_coefficient
        * (
            x - left_boundary
        )
    )
```

For the right tail:

```text
1 - CDF(x) =
    (
        1 - CDF(right_boundary)
    )
    * exp(
        -right_tail_coefficient
        * (
            x - right_boundary
        )
    )
```

The theoretical extrapolation produced more accurate simulated moments than the
empirical boundary estimator.

## Tail-Mass Validation

The probability mass outside the valid FFT region is checked through:

```text
tail_mass =
    max(
        CDF(left_boundary),
        1 - CDF(right_boundary)
    )
```

The selected grid is considered adequate when:

```text
tail_mass <= 0.0001
```

The observed tail masses satisfy this condition for all three models.

## Probability Atom in the Minimal Additive Model

The future increment of the Minimal Additive model has a mixed distribution:

- an absolutely continuous component;
- a probability atom at a deterministic location.

The increment characteristic function contains a non-decaying term:

```text
atom_probability
    * exp(
        imaginary_unit
        * frequency
        * atom_location
    )
```

Because this term does not vanish at high frequencies, applying FFT inversion to
the complete characteristic function would be unstable.

The distribution is therefore decomposed into:

```text
probability atom

continuous conditional distribution
```

Simulation proceeds as follows:

```text
draw Bernoulli variable

if atom is selected:
    increment = atom_location

otherwise:
    simulate from the continuous component
    using Lewis-FFT inversion
```

## Analytical Minimal Additive CDF

The Minimal Additive increment CDF is also derived analytically.

It consists of:

- an exponential left branch;
- a discontinuity at the atom;
- an exponential right branch.

This allows exact inverse-transform simulation without FFT.

Three simulation methods are compared:

```text
FFT with empirically estimated tails

FFT with theoretical tails

analytical CDF inversion
```

The analytical method is the fastest.

The FFT method with theoretical tail coefficients is consistently more accurate
than the version using empirical tail estimates.

For one million simulated increments, the analytical method was approximately
1.7 times faster than the FFT alternatives.

## Forward-Start Option Validation

The Monte Carlo engine is validated using a forward-start call with:

```text
first_reset = 6 months

second_reset = 1 year

strike_multiplier = 1
```

The payoff depends jointly on:

- the additive factor at the first reset;
- the independent increment between the first and second reset.

Five million Monte Carlo scenarios are used in the reported comparison.

The resulting prices are approximately:

```text
Additive Bachelier = 2.9819

Generalized Logistic = 2.9261

Minimal Additive = 2.2942
```

For the Minimal Additive model, an analytical benchmark is available:

```text
analytical_price = 2.29685646

Monte_Carlo_price = 2.29424566
```

The close agreement validates the increment simulation and conditional CDF
reconstruction.

## Part 4 — Exotic Option Pricing

The project prices three exotic products:

- call-on-call;
- put-on-put;
- chooser option.

The relevant dates are:

```text
T1 = 6 months

T2 = 1 year
```

The second strike is set equal to the initial one-year forward:

```text
K2 =
    initial_forward_price(0, T2)
```

For the principal compound-option comparison:

```text
K1 = 0.3
```

## Call-on-Call

The call-on-call gives the holder the right to purchase, at the first reset date,
a call option expiring at the second reset date.

Its value depends on the conditional value of the inner call.

The reported prices are:

```text
Additive Bachelier = 5.5329

Generalized Logistic = 5.5338

Minimal Additive = 5.5303
```

## Put-on-Put

The put-on-put gives the holder an option on the value of a future put.

The reported prices are:

```text
Additive Bachelier = 0.0073

Generalized Logistic = 0.0079

Minimal Additive = 0.0053
```

The selected compound strike produces very small put-on-put values.

## Chooser Option

At the first reset date, the chooser holder selects whether the final payoff will
be a call or a put.

The reported prices are:

```text
Additive Bachelier = 10.3800

Generalized Logistic = 10.3800

Minimal Additive = 10.3780
```

The three models produce similar prices for the selected contracts.

## Avoiding Nested Monte Carlo

The exotic pricing problem naturally leads to nested simulation.

For every outer realization at `T1`, an inner conditional option value at `T2`
must be calculated.

A direct nested Monte Carlo would have computational complexity approximately
equal to:

```text
number_of_outer_paths
    * number_of_inner_paths
```

Two more efficient methods are implemented.

## Stochastic-Mesh Method

The first method:

1. constructs a grid of first-period simulated factors;
2. calculates conditional call and put values on the grid;
3. interpolates the conditional pricing functions;
4. evaluates the interpolants on the outer paths.

The grid combines:

- empirical quantiles;
- central equally spaced points;
- extreme simulated values.

PCHIP interpolation is used because it preserves monotonicity and avoids
artificial oscillations.

## No-Grid Method

The second method avoids interpolation completely.

The second-period increments are sorted once.

Conditional call and put expectations are then computed from:

- binary search;
- cumulative counts;
- cumulative sums.

This produces the same empirical inner expectations without a grid or nested
loops.

The method has approximate complexity:

```text
O(
    number_of_paths
    * log(number_of_paths)
)
```

rather than quadratic nested-simulation complexity.

## No-Grid Performance

The no-grid method is approximately 19 to 24 times faster than the
stochastic-mesh implementation.

Representative speedups are:

```text
Additive Bachelier = 18.97 times

Generalized Logistic = 19.82 times

Minimal Additive = 24.47 times
```

The maximum pricing discrepancy between the two approaches is below:

```text
0.0000001
```

The no-grid method is therefore used for the final exotic prices.

## Inner Put-Call Parity Check

The conditional call and put estimates satisfy empirical put-call parity.

The maximum observed absolute error is approximately:

```text
9.3e-9
```

This confirms the internal consistency of the no-grid implementation.

## Monte Carlo Confidence Intervals

Monte Carlo standard errors and 95% confidence intervals are reported for all
products and models.

The uncertainty is measured both:

- relative to the underlying forward;
- relative to the exotic option price.

For the market-maker quotation test, the forward-based normalization is used.

All 95% confidence-interval half-widths are below 5 basis points of the reference
forward.

The Monte Carlo estimators are therefore sufficiently accurate for a total
bid-ask quotation of 10 basis points of the forward.

## Bid-Ask Quotation

The total exotic bid-ask spread is interpreted as:

```text
10 basis points
    of the reference forward
```

The half-spread is:

```text
5 basis points
    of the reference forward
```

Bid and ask prices are calculated as:

```text
bid =
    max(
        0,
        mid_price
        - half_spread
    )

ask =
    mid_price
    + half_spread
```

The bid is floored at zero because option prices cannot be negative.

## Part 5 — Minimal Additive Analytical Benchmarks

The Minimal Additive model admits elementary formulas for selected exotic
configurations.

Analytical expressions are derived for:

- call-on-call;
- put-on-put;
- chooser option.

The derivations exploit:

- the asymmetric Laplace density of the first increment;
- the probability atom of the second increment;
- exponential continuous tails;
- piecewise-linear option payoffs;
- strategically selected compound strikes.

## Analytical vs Monte Carlo Prices

For the tractable strike configurations, the comparison is:

```text
Call-on-Call:
    Monte Carlo = 4.0264
    Analytical = 4.0292

Put-on-Put:
    Monte Carlo = 0.4903
    Analytical = 0.5088

Chooser:
    Monte Carlo = 10.3788
    Analytical = 10.3807
```

The analytical formulas provide an independent validation of the Monte Carlo
engine.

The put-on-put exhibits the largest discrepancy because it is more sensitive to
interpolation and tail treatment in the numerical conditional-pricing step.

## Part 6 — Risk Management

The project studies the hedging of a short exotic portfolio.

The market maker is assumed to have sold one unit of each product:

```text
short one Call-on-Call

short one Put-on-Put

short one Chooser
```

The risk-management engine uses the Additive Bachelier model because it provides
the lowest global calibration RMSE and a strong fit around the relevant
maturities.

## Hedging Instruments

The allowed hedge instruments are:

- plain vanilla options;
- the front WTI future.

The selected instruments are:

```text
one ATM call close to the 6-month reset

one ATM call close to the 12-month reset

the first future expiring after the valuation date
```

The option strikes are selected as the quoted call strikes closest to the
corresponding market-implied forwards.

Calls are preferred because:

- more call strikes are available near the relevant maturities;
- call-side calibration errors are lower than put-side errors;
- ATM calls provide strong Vega exposure.

## Hedged Risk Factors

The hedge targets three first-order risk factors:

```text
12-month Vega

6-month Vega

forward Delta
```

The hedge is constructed sequentially from the longest volatility exposure to
the shortest.

## Finite-Difference Greeks

Delta and Vega are calculated using central finite differences.

The absolute bumps are:

```text
forward_bump = 0.01

normal_volatility_bump = 0.01
```

Delta is estimated as:

```text
Delta =
    (
        value_with_up_forward
        - value_with_down_forward
    )
    / (
        2 * forward_bump
    )
```

Vega is estimated separately for the 6-month and 12-month volatility pillars.

Common random numbers are reused in the up and down scenarios to reduce Monte
Carlo noise.

## Contractual Strike Treatment

When the forward is bumped, the exotic strike is held fixed.

This is essential because the strike is determined at trade inception and does
not move with subsequent market conditions.

Recalculating the strike after every bump would materially understate the
product's Delta.

## Cascade Hedging Rule

The hedge is constructed in three stages.

### Stage 1 — Hedge 12-Month Vega

The 12-month ATM call is used to neutralize the long-maturity Vega.

```text
quantity_12M_call =
    -portfolio_12M_vega
    / call_12M_vega
```

### Stage 2 — Hedge 6-Month Vega

The 6-month ATM call hedges the remaining short-maturity Vega after accounting
for the 6-month Vega generated by the 12-month call.

```text
quantity_6M_call =
    -remaining_6M_vega
    / call_6M_vega
```

### Stage 3 — Hedge Delta

The front future hedges the residual Delta after the two vanilla options have
been added.

```text
quantity_front_future =
    -residual_portfolio_delta
```

## Integer Hedge Quantities

Hedge positions are restricted to integer quantities.

The continuous hedge ratios are truncated toward zero.

```text
positive_quantity:
    use floor

negative_quantity:
    use ceiling
```

The rule is deliberately conservative:

- it avoids over-hedging;
- it limits transaction costs;
- it leaves a small residual first-order exposure.

The truncation is applied sequentially within the cascade.

## Selected Hedge

The resulting hedge quantities are:

```text
6-month ATM calls = 1

12-month ATM calls = 1

front futures = 0
```

After adding the two calls, the remaining Delta is sufficiently small to be
truncated to zero.

## Hedging Costs

The project assumes full bid-ask spreads of:

```text
front future = 1 basis point

plain vanilla option = 4 basis points
```

Because the hedge is initiated from mid prices, the one-way transaction cost is
one half of the full spread.

The total initial hedging cost is:

```text
0.0004
```

## Initial Cash Flow

The initial cash flow includes:

- premiums received from selling the three exotics;
- premiums paid for the hedging options;
- transaction costs.

The front future has no initial premium.

The resulting initial net cash flow is:

```text
6.2996
```

This cash flow is not included in the reported mark-to-market P&L comparison.

## Hedge Backtest

The hedge is evaluated on:

```text
9 June 2020

16 June 2020
```

On each date:

1. synthetic forwards and discount factors are recalibrated;
2. the volatility surface is rebuilt;
3. the exotic products are repriced;
4. the vanilla hedge options are marked using market prices;
5. the unhedged and hedged P&L are calculated.

## Backtest Results

The reported results are:

```text
9 June 2020:

unhedged P&L = -1.4010

hedged P&L = 0.7587

hedge improvement = 2.1596
```

```text
16 June 2020:

unhedged P&L = -2.5577

hedged P&L = 0.1019

hedge improvement = 2.6596
```

The first-order hedge materially reduces the losses of the short exotic
portfolio on both test dates.

## Sources of Residual Hedging Error

The hedge is not expected to be perfect because it does not explicitly neutralize:

- Gamma;
- Vanna;
- Volga;
- higher-order cross sensitivities;
- model recalibration risk;
- forward-curve basis risk;
- volatility-surface shape risk.

Additional residual errors arise from:

- integer hedge quantities;
- using a front future to hedge a longer-dated forward exposure;
- marking hedge instruments to market while valuing exotics with a model;
- changes in the availability and liquidity of quoted strikes.

## Key Results

The principal numerical findings are:

```text
Best global calibration:
    Additive Bachelier

AB global RMSE:
    0.1007

GL global RMSE:
    0.1052

MA global RMSE:
    0.1336
```

```text
No-grid exotic pricing speedup:
    approximately 20 times
```

```text
Maximum difference between
mesh and no-grid prices:
    below 1e-7
```

```text
MA analytical CDF simulation:
    approximately 1.7 times faster
    than FFT-based simulation
```

```text
Monte Carlo uncertainty:
    compatible with a 10-bps
    forward-based bid-ask quote
```

```text
Hedged P&L:
    substantially more stable
    than unhedged P&L
    on both backtesting dates
```

## Numerical Validation

The implementation includes the following validation checks.

### Synthetic Forward Regression

```text
R-squared close to 1
```

for valid maturities.

### Put-Call Parity

```text
call_price
    - put_price
    approximately equals
discount_factor
    * (
        forward
        - strike
    )
```

### Discount Factor Validity

```text
discount_factor > 0
```

### Characteristic Function at Zero

```text
characteristic_function(0) = 1
```

### Increment Consistency

```text
increment_CF(s, t, u)
    * marginal_CF(s, u)
    approximately equals
marginal_CF(t, u)
```

### CDF Monotonicity

```text
CDF(x_i_plus_1)
    greater than or equal to
CDF(x_i)
```

### Probability Bounds

```text
0 <= CDF <= 1
```

### Tail-Mass Condition

```text
tail_mass <= 1e-4
```

### Minimal Additive Atom

```text
atom_probability
    + continuous_probability
    approximately equals
1
```

### Simulated Moments

The first four simulated moments are compared across:

- empirical-tail FFT;
- theoretical-tail FFT;
- analytical inversion.

### Forward-Start Benchmark

The Minimal Additive Monte Carlo price is compared with its analytical value.

### Inner Put-Call Parity

Conditional no-grid call and put estimates satisfy empirical parity.

### Monte Carlo Confidence Intervals

Confidence intervals are checked against the quotation tolerance.

### Hedge Residuals

Post-hedge Delta and Vega exposures are calculated before and after integer
rounding.

## Computational Design

The MATLAB implementation emphasizes:

- modular model characteristic functions;
- reusable Fourier-pricing functions;
- vectorized Monte Carlo simulation;
- common random numbers;
- precomputed FFT grids;
- theoretical tail extrapolation;
- sorted-increment conditional pricing;
- separation between calibration, simulation, pricing, and risk management.

## Suggested Repository Structure

```text
wti-linear-additive-models/
|
|-- README.md
|
|-- src/
|   |
|   |-- calibration/
|   |   |-- calibrateSyntheticForwards.m
|   |   |-- extractATMVolatility.m
|   |   |-- calibrateAdditiveBachelier.m
|   |   |-- calibrateMinimalAdditive.m
|   |   |-- calibrateGeneralizedLogistic.m
|   |   `-- buildOTMSurface.m
|   |
|   |-- models/
|   |   |-- additiveBachelierLogCF.m
|   |   |-- minimalAdditiveCF.m
|   |   |-- generalizedLogisticLogCF.m
|   |   |-- generalizedLogisticPDF.m
|   |   |-- minimalAdditivePDF.m
|   |   `-- complexLogGamma.m
|   |
|   |-- fourier/
|   |   |-- priceEuropeanLewisFFT.m
|   |   |-- invertCDFLewisFFT.m
|   |   |-- restrictMonotoneCDF.m
|   |   |-- extrapolateCDFTails.m
|   |   `-- invertNumericalCDF.m
|   |
|   |-- simulation/
|   |   |-- simulateAdditiveIncrement.m
|   |   |-- simulateMinimalAdditiveIncrement.m
|   |   |-- simulateMinimalAdditiveAnalytical.m
|   |   `-- simulateAdditivePath.m
|   |
|   |-- pricing/
|   |   |-- priceForwardStartOption.m
|   |   |-- priceCompoundOptionsNoGrid.m
|   |   |-- priceCompoundOptionsMesh.m
|   |   |-- priceChooserOption.m
|   |   `-- conditionalVanillaValues.m
|   |
|   |-- analytical/
|   |   |-- minimalAdditiveIncrementCDF.m
|   |   |-- minimalAdditiveForwardStart.m
|   |   |-- minimalAdditiveCallOnCall.m
|   |   |-- minimalAdditivePutOnPut.m
|   |   `-- minimalAdditiveChooser.m
|   |
|   |-- risk_management/
|   |   |-- computeFiniteDifferenceGreeks.m
|   |   |-- selectHedgingInstruments.m
|   |   |-- buildCascadeHedge.m
|   |   |-- computeHedgingCost.m
|   |   `-- backtestHedge.m
|   |
|   `-- utilities/
|       |-- bachelierPrice.m
|       |-- bachelierImpliedVolatility.m
|       |-- calculateConfidenceInterval.m
|       |-- calculateRMSE.m
|       `-- yearFractionACT365.m
|
|-- scripts/
|   |-- runForwardCalibration.m
|   |-- runModelCalibration.m
|   |-- runSimulationValidation.m
|   |-- runExoticPricing.m
|   |-- runAnalyticalBenchmarks.m
|   `-- runRiskManagement.m
|
|-- data/
|   `-- README.md
|
|-- results/
|   |-- calibration/
|   |-- simulation/
|   |-- exotic_pricing/
|   |-- hedging/
|   `-- figures/
|
`-- report/
    `-- final_project_report.pdf
```

The exact names can be adapted to the structure of the original MATLAB library.

## Requirements

The project was developed in MATLAB.

Representative requirements include:

```text
MATLAB

Optimization Toolbox

Statistics and Machine Learning Toolbox
```

Some functions may be implemented without additional toolboxes depending on the
numerical routines used.

## Running the Project

A recommended workflow is:

```text
1. Add the src directory and its subfolders to the MATLAB path.

2. Configure the market-data directory.

3. Run the synthetic-forward calibration.

4. Run the three model calibrations.

5. Generate the increment CDFs and simulation tests.

6. Price the forward-start and exotic products.

7. Run the Minimal Additive analytical benchmarks.

8. Run the risk-management and hedging backtest.
```

Representative MATLAB commands are:

```matlab
addpath(genpath("src"));

run("scripts/runForwardCalibration.m");
run("scripts/runModelCalibration.m");
run("scripts/runSimulationValidation.m");
run("scripts/runExoticPricing.m");
run("scripts/runAnalyticalBenchmarks.m");
run("scripts/runRiskManagement.m");
```

## Main Outputs

The project produces:

- market-implied discount factors;
- WTI forward curve;
- continuously compounded zero rates;
- deterministic absolute dividends;
- ATM normal-volatility term structure;
- calibrated AB, MA, and GL parameters;
- calibration errors by maturity and moneyness;
- model PDFs and tail comparisons;
- numerical increment CDFs;
- simulated additive paths;
- analytical and simulated Minimal Additive distributions;
- forward-start option prices;
- call-on-call prices;
- put-on-put prices;
- chooser-option prices;
- Monte Carlo standard errors;
- bid and ask quotations;
- analytical Minimal Additive exotic prices;
- Delta and bucketed Vega sensitivities;
- hedge quantities;
- transaction costs;
- unhedged and hedged backtesting P&L.

## Technologies and Topics

- MATLAB
- WTI crude oil options
- Commodity derivatives
- Additive Bachelier model
- Minimal Additive model
- Generalized Logistic model
- Synthetic forwards
- Normal implied volatility
- Characteristic functions
- Lewis Fourier inversion
- Fast Fourier Transform
- Monte Carlo simulation
- Mixed distributions
- Probability atoms
- Compound options
- Chooser options
- Numerical calibration
- Delta–Vega hedging
- Quantitative risk management

## Data Availability

The original project uses course-provided WTI option data.

The public repository should not contain restricted or proprietary market data
unless redistribution is explicitly permitted.

When the original files cannot be published, the `data` directory should contain
a description of:

- expected folder structure;
- option file naming convention;
- required columns;
- date format;
- strike and price units;
- call and put identifiers;
- futures-expiry mapping;
- missing-data treatment.

Synthetic or anonymized sample data may be included to demonstrate the expected
input format.

## Reproducibility

For reproducible Monte Carlo results, the project should report:

```text
random seed

number of simulations

FFT grid size

FFT range

complex shift rule

tail extrapolation method

optimizer initial values

optimizer tolerances

finite-difference bumps
```

The reported results use fixed random seeds and common random numbers where
appropriate.

## Academic Context

This project was developed as the final project of the Financial Engineering
course at Politecnico di Milano.

Authors:

```text
Anna Belli

Gabriele Bulian

Marco Gretter
```

The repository presents the MATLAB library and the complete quantitative workflow
developed for the project:

```text
market calibration
    -> model calibration
    -> Fourier simulation
    -> exotic pricing
    -> analytical validation
    -> hedging and backtesting
```

## References

The methodology builds on academic research concerning:

- fast Monte Carlo simulation for additive processes;
- synthetic forwards and market-implied funding curves;
- the Additive Bachelier model for oil options;
- additive logistic processes;
- Fourier simulation of Lévy-driven and additive models.

Full bibliographic references are provided in the project report.
