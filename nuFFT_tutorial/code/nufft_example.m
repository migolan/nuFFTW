% NUFFT_EXAMPLE - Usage example of the nufft toolbox
%
%   The nuFFT implementation here is based on the paper
%   "Rapid Gridding Reconstruction With a Minimal Oversampling Ratio",
%   P.J. Beatty, D.G. Nishimura and J.M. Pauly, IEEE Transactions on
%   Medical Imaging, Vol. 24, No. 6, June 2005.
%
%   This implementation assumes 2-dimensional data with equal gridding
%   parameters in both axes.
%
%   The interpolation kernel used for gridding is the Kaiser-Bessel kernel.

% Michal Zarrouk, June 2013.


%% load data and initialize parameters

addpath ../../utils
load ../data/noncartesian_phantom.mat
% This data file contains the following variables:
% kcoord - k-space coordinates of non-uniform 2d samples
% dcf    - density compensation factors
% kdata  - k-space data samples
% imagedirect - exact non-uniform FT of the k-space samples, on 128x128
%   image, computed the slow way (computation takes about a minute using
%   one-threaded C++ code on my MacBook Pro)

sqrtdcf = sqrt(dcf);
N = 128; % output image size
maxerr = 1e-3; % maximum aliasing error allowed, as defined in Beatty et al. formula (3).
alpha = 2; % grid oversampling ratio

% This stage might take a while, because of the convolution matrix
% computation. I suggest that once you've computed it, save it for future
% reconstructions.
tic
nufft_st = init_nufft(kcoord, sqrtdcf, N, maxerr, alpha);
toc

%%
% adjoint nufft
kdata_denscomp = kdata.*sqrtdcf;
imagenufft = adjoint_nufft(kdata_denscomp, nufft_st);
e = get_image_errors(imagenudft,imagenufft);
max(e(:))
figure
imshowz(e)

% forward nufft
kdata_denscomp2 = forward_nufft(imagenufft, nufft_st);

% check data consistency:
kdata_denscomp'*kdata_denscomp2-imagenufft(:)'*imagenufft(:)

%%
% lets go backward again:
imagenufft2 = adjoint_nufft(kdata_denscomp2, nufft_st);

% and check data consistency in the other direction:
imagenufft(:)'*imagenufft2(:) - kdata_denscomp2'*kdata_denscomp2

