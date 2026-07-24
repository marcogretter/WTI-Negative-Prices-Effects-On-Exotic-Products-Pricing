function price = price_Bachelier_Lewis_FFT(x_target,t, B0,a, sigma_t, eta, kappa, M, value, type, model)
% European Call pricing via Lewis formula using FFT
%
% Inputs:
%   x_target - Vector of Bachelier moneyness target points.
%              In the original pricing formula: x = K - F0.
%              In the normalized framework: x_target is y.
%
%   t        - Time to maturity. For normalized pricing use t = 1.
%   B0       - Discount factor. For normalized pricing use B0 = 1.
%   a        - Lewis shift, must lie inside the analyticity strip.
%   sigma_t  - Model volatility level. For normalized pricing use sigma_t = 1.
%   eta      - Skew parameter.
%   kappa    - Vol-of-vol parameter.
%   M        - Grid size exponent, N = 2^M.
%   value    - Numerical value of dz or x1.
%   type     - Either 'dz' or 'x1'.
%
% Output:
%   price    - Call prices corresponding to x_target.
%
% Formula implemented:
%
% C(x,t) = B0 * [ R_a(x)
%                 + exp(a*x)/(2*pi) * int phi(xi+i*a)
%                   exp(-i*xi*x)/(i*xi-a)^2 dxi ]
%
% where:
%   R_a(x) = 0    if a < 0
%   R_a(x) = -x   if a > 0

% 1. Parameters & Relations
N = 2^M;

if nargin < 11
    model = 'AB';
end

% Force column vector
x_target = x_target(:);

% Grid Logic based on selection
if strcmpi(type, 'dz')
    % User selected dz
    dz = value;
    dx = (2*pi) / (N * dz);
    x1 = -dx * (N-1) / 2;

elseif strcmpi(type, 'x1')
    % User selected x1, i.e. left endpoint of Fourier grid
    x1 = value;
    dx = -2 * x1 / (N-1);
    dz = (2*pi) / (N * dx);

else
    error('Invalid type. Choose either ''dz'' or ''x1''');
end

% Sort target points
[x_target_sorted, idx] = sort(x_target);

% Symmetry of the grids
x1 = -dx * (N-1) / 2;
xN = -x1;
x_grid = (x1 : dx : xN).';

z1 = -dz * (N-1) / 2;
zN = -z1;
z_grid = (z1 : dz : zN).';

% Safety: force exact length N if colon creates N+1 because of floating point
x_grid = x_grid(1:N);
z_grid = z_grid(1:N);

% 2. Definition of Integrand
u = x_grid + 1i*a;

switch upper(model)
    case 'AB'
        cf_values = additiveBachelierCf(u, t, sigma_t, eta, kappa);

    case 'GL'
        % Here eta is used as alpha_GL and kappa as beta_GL.
        alphaGL = eta;
        betaGL  = kappa;

        cf_values = generalizedLogisticCf(u, alphaGL, betaGL);

    otherwise
        error('Unknown model. Use ''AB'' or ''GL''.');
end


f_x = (1/(2*pi)) .* cf_values ./ ((1i*x_grid - a).^2);

% 3. FFT
j_minus_1 = (0:N-1).';

input_fft = f_x .* exp(-1i * z1 * dx * j_minus_1);

% Fast Fourier Transform
Y = fft(input_fft);

% Integral Reconstruction via Prefactor
integral_values = dx .* exp(-1i * x1 .* z_grid) .* Y;

% 4. Interpolation to target moneyness
integral_interp = interp1(z_grid, real(integral_values), x_target_sorted, 'spline', NaN);

% 5. Residual term R_a(x)
if a > 0
    R = -x_target_sorted;
elseif a < 0
    R = zeros(size(x_target_sorted));
else
    error('Use a non-zero shift a.');
end

% 6. Full Lewis price
prices_sorted = B0 .* (R + exp(a .* x_target_sorted) .* integral_interp);

% Prices in the original order
price = NaN(size(x_target_sorted));
price(idx) = prices_sorted;

end