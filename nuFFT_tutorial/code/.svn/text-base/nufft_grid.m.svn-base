function [gridded_data, convmat] = nufft_grid (kcoord, sqrtdcf, kdata, N, maxerr, alpha, W)
% NUFFT_GRID - Grid non-uniformly sampled data onto a uniform grid
%   
%   gridded_data = NUFFT_GRID (kcoord, kdata, N, maxerr, alpha, W)
%
%   Arguments:
%   kcoord - k-sapce sample coordinates, normalized to region [-0.5,0.5].
%       Array size is [Mx2], where M is the number of nonuniform samples.
%   sqrtdcf - square root of density compensation factors.
%   kdata - density compensated non-uniformly sampled k-space data, array
%       length is M.
%   N (optional) - Image size in pixels. Default value is 128.
%   maxerr (optional) - Maximum aliasing error allowed, as defined in
%       Beatty et al. Default value is 1e-3.
%   alpha (optional) - Oversampling ratio. Default value is 2.
%   W (optional) - Interpolation kernel width. Default value is ~3.84,
%       empirically found to match alpha=2 and maxerr=1e-3, according to
%       Beatty's paper.
%   gridded_data - data gridded to a uniform ovesampled grid
%
%   This function assumes 2-dimensional data with equal gridding parameters
%   in both axes.
%   The interpolation kernel used is the Kaiser-Bessel, as defined in
%   Beatty et al.

% Michal Zarrouk, June 2013.


[convmat, gridded_data] = compute_gridding_matrix(kcoord, sqrtdcf, N, maxerr, alpha, W, 0, kdata);
