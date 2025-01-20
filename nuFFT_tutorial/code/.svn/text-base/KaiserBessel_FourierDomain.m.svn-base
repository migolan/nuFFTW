function KB = KaiserBessel_FourierDomain (k, W, G, beta)
% KaiserBessel_FourierDomain - Kaiser-Bessel function in the Fourier domain
%   as presented by Beatty et al. in formula (4).
% 
%   KB = KaiserBessel_FourierDomain (k, W, G, beta)
%
%   Arguments:
%   k - spatial frequency (units 1/pixel)
%   W - Kaiser-Bessel window width in pixel units
%   G - oversampled image size
%   beta - Kaiser-Bessel window shape factor
%   KB - Kaiser-Bessel function in the Fourier domain

% Michal Zarrouk, June 2013.


KB = G/W*besseli(0,beta*sqrt(1-(2*k*G/W).^2));

% k-space kernel bounds, in k-space units (1/pixel)
kerrad = W/2/G;

% find k-space locations that are outside of the kernel bounds and zero
% them out
KB(abs(k) > kerrad) = 0;
