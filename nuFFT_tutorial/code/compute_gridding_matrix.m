function [convmat, gridded_data] = compute_gridding_matrix(kcoord, sqrtdcf, N, maxerr, alpha, W, phasemod, kdata)
% COMPUTE_GRIDDING_MATRIX - Computes the gridding convolution matrix
%   from a nonuniform sampling pattern onto a uniform oversampled grid.
%   
%   [convmat, gridded_data] = COMPUTE_GRIDDING_MATRIX(kcoord, N, maxerr,...
%       alpha, W, phasemod, kdata)
%
%   Arguments:
%   kcoord - k-sapce sample coordinates, normalized to region [-0.5,0.5].
%       Array size is [Mx2], where M is the number of nonuniform samples.
%   sqrtdcf - square root of density compensation factors.
%   N (optional) - Image size in pixels in each dimension. Default value is 128.
%   maxerr (optional) - Maximum aliasing error allowed, as defined in
%       Beatty et al. Default value is 1e-3.
%   alpha (optional) - Oversampling ratio. Default value is 2.
%   W (optional) - Interpolation kernel width. Default value is 3.84,
%       empirically found to match alpha=2 and maxerr=1e-3, according to
%       Beatty's paper.
%   phasemod - (optional) flag stating whether to apply phase modulation to
%       avoid fftshifting after the fft (default is 1).
%   kdata (optional) - sampled k-space data from the nonuniform sample
%       locations kcoord, array length is M.
%   convmat - gridding (convolution) matrix, from nonuniform to uniform
%       grid (for adjoint nufft). Array size is [G^2xM], where G is the
%       oversampled grid size.
%   gridded_data - (optional) If sampled data was inserted (as kdata), this
%       will be the gridded data on a uniform ovesampled grid.
%
%   This function assumes 2-dimensional data with equal gridding parameters
%   in both axes.
%   The interpolation kernel used is the Kaiser-Bessel, as defined in
%   Beatty et al.

% Michal Zarrouk, June 2013.


% Set default parameter values
if nargin < 7 || isempty(phasemod)
    phasemod = 1;
end
if nargin < 6 || isempty(W)
    W = 3.84326171875 * ones(1,2); % (empirically found to match alpha=2 and maxerr=1e-3, according to Beatty's paper)
end
if nargin < 5 || isempty(alpha)
    alpha = 2;
end
if nargin < 4 || isempty(maxerr)
    maxerr = 1e-3;
end
if nargin < 3 || isempty(N)
    N = 128 * ones(1,2);
end

% number of pixels in each dimension of the oversampled image
G = OversampledImageSize(alpha,N);

%% presample the interpolation kernel
% we need to sample the kernel densly enough to not ruin the nuFFT - this
% is the minimial kernel presampling ratio to satisfy the accuracy we want,
% according to Beatty's paper, formulas 7 and 8.

% Kaiser-Bessel kernel shape parameter, best suited for maxerr, according
% to Beatty's paper
beta = KaiserBesselOptimalBeta (W, alpha);

% k-space kernel radius, in k-space units (1/pixel)
kerrad = W/2./G;

% presampled kernel interpolation method: 'NEAREST_NEIGHBOR' or 'LINEAR'
kernel_interp_method = 'LINEAR';

% sampling ratio of the gridding interpolation kernel
S = KernelSamplingRatio (maxerr, alpha, kernel_interp_method);

% increment of presampled kernel samples, in k-space units (1/pixel)
kerinc = 1/S./G;

for dim = 1:length(N)
    % k-space locations of presampled kernel
    k_kb{dim} = (-kerrad(dim)-kerinc(dim)):kerinc(dim):(kerrad(dim)+kerinc(dim));

    % presampled Kaiser-Bessel kernel, formula 4 in Beatty's paper
    KB{dim} = KaiserBessel_FourierDomain (k_kb{dim}, W(dim), G(dim), beta(dim));
end
%b = pi*(2-1/alpha);
%KB = sinh(b*sqrt((W/2)^2-(G*k_kb).^2))./sqrt((W/2)^2-(G*k_kb).^2)/pi;
%KB(a) = sin(b*sqrt((G*k_kb(a)).^2-(W/2)^2))./sqrt((G*k_kb(a)).^2-(W/2)^2)/pi;

%% compute values
% uniform k-space grid coordinates, oversampled by alpha
ksx = -0.5:(1/G(1)):0.5;
ksy = -0.5:(1/G(2)):0.5;
[ksy,ksx] = meshgrid(ksy(1:G(1)),ksx(1:G(2)));

M = size(kcoord,1); % number of nonuniform k-space samples
nnz = 0; % number of nonzeros in the matrix
nnz_est = ceil(prod(W)*M); % estimate of nnz
% pre-allocate sparse gridding (convolution) matrix
row_idcs = zeros(nnz_est,1);
col_idcs = zeros(nnz_est,1);
values   = zeros(nnz_est,1);

% If we have kdata to grid:
if nargin > 6 && nargout > 1
    gridded_data = zeros(size(G));
end

% phase modulation to avoid fftshifting after the fft
if phasemod
    if mod(G,2) == 0
        % this creates a checkerboard of +-1
        [ax,ay] = meshgrid(1:G(1),1:G(2));
        mask = (-1) .^ (ax+ay);
    else
        % in this case the phase modulation is complex and thus requires double
        % storage for the sparse convolution matrix, and complex multiplication
        % when gridding. usually we would prefer to just use an even G.
        ax = -(0:(G(1)-1))*pi/G(1);
        ax(2:2:end) = ax(2:2:end)+pi;
        ay = -(0:(G(2)-1))*pi/G(2);
        ay(2:2:end) = ay(2:2:end)+pi;
        [ax,ay] = meshgrid(ax,ay);
        mask = exp(1i*(ax+ay));
    end
    mask = mask(:);
end

% for every noncartesian sample
for i = 1:M
    % find samples within kernel width
    nxy = abs(ksx-kcoord(i,1))<=kerrad(1) & abs(ksy-kcoord(i,2))<=kerrad(2);

    % Interpolation of kernel values at grid locations from presampled
    % kernel. Interpolation is performed independently in each axis (i.e.
    % square kernel)
    switch kernel_interp_method
        case 'LINEAR' % linear interpolation of presampled kernel
            conv_x = interp1(k_kb{1}, KB{1}, ksx(nxy)-kcoord(i,1));
            conv_y = interp1(k_kb{2}, KB{2}, ksy(nxy)-kcoord(i,2));
        case 'NEAREST_NEIGHBOR' % nearest-neighbor interpolation of presampled kernel
            % TBD, but you should really just use linear interpolation
            % because it requires less presampling.
    end    
    
    idcs = find(nxy);
    nidcs = length(idcs);
    new_loc = nnz+(1:nidcs);
    row_idcs(new_loc) = idcs;
    col_idcs(new_loc) = i;
    values(new_loc) = conv_x .* conv_y * sqrtdcf(i);
    if phasemod
        values(new_loc) = values(new_loc) .* mask(nxy);
    end
    nnz = nnz+nidcs;
    
    % If we have kdata to grid
    if nargin > 6 && nargout > 1
        gridded_data(nxy) = gridded_data(nxy) + values(new_loc) * kdata(i);
    end
end
convmat = sparse(row_idcs(1:nnz), col_idcs(1:nnz), values(1:nnz), G^2, M);

