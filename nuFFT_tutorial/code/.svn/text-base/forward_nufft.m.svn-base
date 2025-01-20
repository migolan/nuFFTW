function kdata = forward_nufft(im, nufft_st)
% FORWARD_NUFFT - Forward nufft operator
%
%   kdata = FORWARD_NUFFT(im, nufft_st)
% 
%   Arguments:
%   im - Image [NxN].
%   nufft_st - data structure defining the non-uniform Fourier transform
%       and it's adjoint, between k-space data sampled on the non-uniform
%       sampling locations kcoord, and a uniform grid in image space.
%   kdata - k-space data resulting from non-uniform Fourier transform of
%       the image.
%
%   This function assumes 2-dimensional data with equal gridding parameters
%   in both axes.
%   See nufft_example.m for a usage example.

% Michal Zarrouk, June 2013.


% deapodize
% our deapodization matrix is actually real and symmetric so we don't really
% need to to conj(nufft_st.deapmat)
im_deap = im./conj(nufft_st.deapmat)/nufft_st.G^2;

% zero-pad for oversampled image
im_oversampled = zeros(nufft_st.G);
im_oversampled(nufft_st.nstart+(1:nufft_st.N),nufft_st.nstart+(1:nufft_st.N)) = im_deap;

% Fourier transform
gridded_data = fft2(im_oversampled);

% convolve (grid) uniform data onto nonuniform k-space locations
gridded_data = gridded_data(:);
kdata = nufft_st.convmat' * gridded_data;

