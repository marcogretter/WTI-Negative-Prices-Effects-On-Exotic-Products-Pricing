function pdf = pdf_Bachelier_FFT(y_target, eta, kappa, M, value, type)
% PDF of normalized Additive Bachelier variable Y via Fourier inversion
%
% Computes:
%   p(y) = (1/(2*pi)) int exp(-i*z*y) phi_Y(z) dz
%
% where:
%   phi_Y(z) = additiveBachelierCf(z, 1, 1, eta, kappa)
%
% Inputs:
%   y_target - target points where the pdf is evaluated
%   eta      - calibrated / chosen eta
%   kappa    - calibrated / chosen kappa
%   M        - grid exponent, N = 2^M
%   value    - numerical value of dz or x1
%   type     - 'dz' or 'x1'
%
% Output:
%   pdf      - pdf values at y_target

N = 2^M;

% Force column vector
y_target = y_target(:);

% 1. Grid Logic based on selection
if strcmpi(type, 'dz')
    dz = value;
    dy = (2*pi) / (N * dz);
    y1 = -dy * (N-1) / 2;

elseif strcmpi(type, 'x1')
    y1 = value;
    dy = -2 * y1 / (N-1);
    dz = (2*pi) / (N * dy);

else
    error('Invalid type. Choose either ''dz'' or ''x1''');
end

% Sort target points
[y_target_sorted, idx] = sort(y_target);

% Symmetric grids
y1 = -dy * (N-1) / 2;
yN = -y1;
y_grid = (y1 : dy : yN).';

z1 = -dz * (N-1) / 2;
zN = -z1;
z_grid = (z1 : dz : zN).';

% Safety in case colon gives N+1 points
y_grid = y_grid(1:N);
z_grid = z_grid(1:N);

% 2. Fourier integrand for PDF
% No Lewis denominator, no shift a.
f_z = (1/(2*pi)) .* additiveBachelierCf(z_grid, 1, 1, eta, kappa);

% 3. FFT reconstruction
j_minus_1 = (0:N-1).';

input_fft = f_z .* exp(-1i * y1 * dz * j_minus_1);

Y = fft(input_fft);

density_grid_complex = dz .* exp(-1i * z1 .* y_grid) .* Y;

density_grid = real(density_grid_complex);

% 4. Interpolate to target points
pdf_sorted = interp1(y_grid, density_grid, y_target_sorted, 'spline', NaN);

% Restore original order
pdf = NaN(size(y_target_sorted));
pdf(idx) = pdf_sorted;

% Numerical cleanup: tiny negative values due to FFT oscillations
pdf(pdf < 0 & pdf > -1e-10) = 0;

end