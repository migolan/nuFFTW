% NUFFT_TUTORIAL - a hands-on nuFFT tutorial
%
%   The purpose of this tutorial is to show the effects of each stage in
%   the nufft (gridding, density compensation, deapodization,
%   oversampling/FOV extraction, fftshifting/+-1 masking, etc.).
%
%   If you just want to learn how to use this nufft toolbox, refer to
%   nufft_example.m.
%
%   The nuFFT implementation here is based on Beatty, Nishimura and Pauly,
%   "Rapid gridding reconstruction with a minimal oversampling ratio".
%
%   It assumes 2-dimensional data with equal gridding parameters in both
%   axes.

% Michal Zarrouk, June 2013.

%% load the data and set nuFFT parameters

addpath ../../utils
load ../data/noncartesian_phantom.mat
% This data file contains the following variables:
% kcoord - k-space coordinates of non-uniform 2d samples
% dcf    - density compensation factors
% kdata  - k-space data samples
% imagedirect - exact non-uniform FT of the k-space samples, on 128x128
%   image, computed the slow way (computation takes about a minute using
%   one-threaded C++ code)

N = 128; % output image size (in pixels)
maxerr = 1e-3; % maximum aliasing error allowed, as defined be in Beatty's paper, formula 3.
alpha = 2; % grid oversampling ratio

% gridding (interpolation) kernel width (empirically found to match alpha
% and maxerr, according to Beatty's paper, see fig. 3)
W = KernelWidth(maxerr, alpha, N);

% let's look at the sample locations
h1 = figure;
plot(kcoord(:,1),kcoord(:,2),'.') % it's a spiral
xlabel('k_x')
ylabel('k_y')

%% convolve the non-uniform data onto a uniform grid and Fourier transform the gridded data
gridded_data = nufft_grid(kcoord, ones(size(kdata)), kdata, N, maxerr, alpha, W);
image_nodeap = ifft2(gridded_data);
imshowz(image_nodeap)

%% oops, we didn't make sure the zero-frequency component is first
% nevermind, let's just fftshift. not centering the data only affects the
% phase, and here we will only look at the magnitude.
image_nodeap = ifftshift(image_nodeap);
imshowz(image_nodeap)

%% We oversampled, so now we have to extract the field of view
G = ceil(alpha*N); % oversampled grid size
% nstart is the index within the oversampled image where the FOV starts.
nstart = ceil((G-N)/2);

image_nodeap = image_nodeap(nstart+(1:N),nstart+(1:N));

subplot(2,2,1)
imshowz(imagedirect)
title('direct non-uniform FT','fontsize',14)
subplot(2,2,2)
imshowz(image_nodeap)
title({'nuFFT without density compensation','or deapodization'},'fontsize',14)
% it's blurred because we didn't compensate for the non-uniform sampling
% density

h2 = figure;
c = floor(N/2);
plot(abs(normalize(imagedirect(c,:))))
hold on
plot(abs(normalize(image_nodeap(c,:))),'r')

%% with density compensation
sqrtdcf = sqrt(dcf);
kdata_denscomp = kdata .* sqrtdcf;
[gridded_data, convmat] = nufft_grid(kcoord, sqrtdcf, kdata_denscomp, N, maxerr, alpha, W);
image_nodeap = ifftshift(ifft2(gridded_data));
image_nodeap = image_nodeap(nstart+(1:N),nstart+(1:N));

figure(h1)
subplot(2,2,3)
imshowz(image_nodeap)
title({'nuFFT with density compensation,','without deapodization'},'fontsize',14)

figure(h2)
plot(abs(normalize(image_nodeap(c,:))),'g')

% looks better, but the interpolation introduces weighting to the data,
% which we haven't compensated for. That's called deapodization.

%% deapodization
% The easiest way to do this is to grid an impulse - 1 in the origin and
% zero elsewhere. The other way is to just calculate the FT of the
% interpolation kernel.

gridded_impulse = nufft_grid([0 0], 1, 1, N, maxerr, alpha, W);
deap = ifftshift(ifft2(gridded_impulse));
deap = deap(nstart+(1:N),nstart+(1:N));
h3 = figure;
imshowz(deap)
title('deapodization function for KB kernel','fontsize',14)

image_deap = image_nodeap ./ deap;

figure(h1)
subplot(2,2,4)
imshowz(image_deap)
title({'nuFFT with density compensation','and deapodization'},'fontsize',14)

figure(h2)
plot(abs(normalize(image_deap(c,:))),'k')
h = legend('direct non-uniform FT','nuFFT, no dens.comp, no deap.','nuFFT, with dens.comp, no deap.','nuFFT, with dens.comp, with deap.');
set(h,'fontsize',14)


%% let's look at the errors
e = get_image_errors(imagedirect,image_deap);
max(e(:))
figure
imshowz(e)
colorbar

%% that was the adjoint nufft, now we go forward

% deapodize
%   It's not intuitive, but going backward actually means dividing by the
% deapodization factors again - because we're not inverting the operations,
% we're just implementing the adjoint, which is just multiplying by the
% transpose conjugate.
%   In our case, the deapodization matrix would be
% D = diag(1./deap), so it's adjoint is D* = diag(1./deap').
im2 = image_deap ./ deap';

% zero-pad for oversampled image
im3 = zeros(G);
im3(nstart+(1:N),nstart+(1:N)) = im2;

% Fourier transform
gridded_data2 = fft2(fftshift(im3));

% we have to compensate for the N factor introduced by the Fourier
% transform
gridded_data2 = gridded_data2/G^2;

% convolve (grid) uniform data onto nonuniform k-space locations
gridded_data2 = gridded_data2(:);
kdata_denscomp2 = convmat' * gridded_data2;

% check consistency of forward and adjoint nufft operators
kdata_denscomp'*kdata_denscomp2 - image_deap(:)'*image_deap(:)

%% nuFFT said.



