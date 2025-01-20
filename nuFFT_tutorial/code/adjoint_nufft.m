function im = adjoint_nufft(kdata, nufft_st)
% ADJOINT_NUFFT - Adjoint nufft operator
%
%   im = ADJOINT_NUFFT(kdata, nufft_st)
% 
%   Arguments:
%   kdata - sampled k-space data (density compensated).
%   nufft_st - data structure defining the non-uniform Fourier transform
%       and it's adjoint, between k-space data sampled on the non-uniform
%       sampling locations kcoord, and a uniform grid in image space.
%   im - Image created by adjoint non-uniform Fourier transform of kdata,
%       array size is [NxN].
%
%   This function assumes 2-dimensional data with equal gridding parameters
%   in both axes.
%   See nufft_example.m for a usage example.

% Michal Zarrouk, June 2013.


kdata = kdata(:);

% grid (convolve) the non-uniform k-space data onto a uniform grid
gridded_data = nufft_st.convmat * kdata;
gridded_data = reshape(gridded_data, nufft_st.G, nufft_st.G);

% a different way to do this would be
% gridded_data = nufft_grid(kcoord, kdata, N, alpha, W, beta, maxerr);

% inverse Fourier transform
image_nodeap = ifft2(gridded_data);

% truncate margins and leave FOV
image_nodeap = image_nodeap(nufft_st.nstart+(1:nufft_st.N),nufft_st.nstart+(1:nufft_st.N));

% deapodize
im = image_nodeap ./ nufft_st.deapmat;

