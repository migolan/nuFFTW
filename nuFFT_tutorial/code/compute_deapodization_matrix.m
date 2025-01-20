function deapmat = compute_deapodization_matrix(N, maxerr, alpha, W, phasemod)
% COMPUTE_DEAPODIZATION_MATRIX - Computes the gridding deapodization matrix
%   from a nonuniform sampling pattern onto a uniform oversampled grid.
%
%   deapmat = COMPUTE_DEAPODIZATION_MATRIX(N, maxerr, alpha, W, phasemod)
%
%   Arguments:
%   N (optional) - Image size in pixels in each dimension. Default value is 128.
%   maxerr (optional) - Maximum aliasing error allowed, as defined in
%       Beatty et al. Default value is 1e-3.
%   alpha (optional) - Oversampling ratio. Default value is 2.
%   W (optional) - Interpolation kernel width. Default value is 3.84,
%       empirically found to match alpha=2 and maxerr=1e-3, according to
%       Beatty's paper.
%   phasemod - (optional) flag stating whether to apply phase modulation to
%       avoid fftshifting after the fft (default is 1).
%   deapmat - Deapodization factors for the image FOV, array size is [NxN].
%
%   This function assumes 2-dimensional data with equal gridding parameters
%   in both axes.
%   The interpolation kernel used is the Kaiser-Bessel, as defined in
%   Beatty's paper.

% Michal Zarrouk, June 2013.


% Set default parameter values
if nargin < 5 || isempty(phasemod)
    phasemod = 1;
end
if nargin < 4 || isempty(W)
    W = 3.84326171875 * ones(1,2); % (empirically found to match alpha=2 and maxerr=1e-3, according to Beatty's paper)
end
if nargin < 3 || isempty(alpha)
    alpha = 2;
end
if nargin < 2 || isempty(maxerr)
    maxerr = 1e-3;
end
if nargin < 1 || isempty(N)
    N = 128 * ones(1,2);
end

%%
% number of pixels in each dimension of the oversampled image
G = OversampledImageSize(alpha,N);

% Kaiser-Bessel kernel shape parameter, best suited for maxerr, according
% to Beatty's paper
beta = KaiserBesselOptimalBeta (W, alpha);

% We can compute the deapodization factors analytically, using the
% expression for our gridding kernel in the image domain (see formula 4 in
% Beatty's paper).
x = -floor(N/2):ceil(N/2-1); % image pixel coordinates
deap = KaiserBessel_ImageDomain (x, W, G, beta); % one-dimensional deapodization factors

% correction for linear interpolation from a presampled gridding kenrel
% optimal sampling ratio for the gridding interpolation kernel
S = KernelSamplingRatio (maxerr, alpha, 'LINEAR');
deap = deap .* sinc(x/S/G).^2;
%SG = ceil(S*G);
%deap1 = deap .* (-1).^(1:N).*sinc(x/SG).^2;

% phase modulation to avoid fftshifting before the fft
if phasemod
    deap = deap .* (-1).^(1:N);
end

deapmat = deap'*deap; % two-dimensional deapodization factors


%   Another way to compute the deapodization factors is to grid an impulse
% (1 in the origin and zero elsewhere).
%   The staright-forward way to grid an impulse would be compute the
% "gridding matrix" of a data point in the origin ([0 0]) and then multiply
% it by a kdata = 1. But this means that for a one-sample kdata (kdata=1),
% the gridding matrix is in fact the gridded data.
%   For some reason, this method produces a slightly bigger error (really
% unnoticeable), and exhibits slightly less consistency between the forward
% and adjoint operators, so I prefer to use this one.
%
% function deapmat = compute_deapodization_matrix(N, alpha, W, maxerr)
% % oversampled grid size
% G = OversampledImageSize(alpha,N);
% % nstart is the index within the oversampled image where the FOV starts.
% nstart = ceil((G-N)/2);
% gridded_impulse = full(compute_gridding_matrix([0 0], N, alpha, W, maxerr));
% gridded_impulse = reshape(gridded_impulse,G,G);
% deapmat = ifft2(gridded_impulse);
% deapmat = deapmat(nstart+(1:N),nstart+(1:N));
%
% This is equivalent to:
% gridded_impulse = nufft_grid([0 0], 1, N, alpha, W, maxerr);
