function kb = KaiserBessel_ImageDomain (x, W, G, beta)
% KaiserBessel_ImageDomain - Kaiser-Bessel function in the image domain
%   as presented by Beatty et al. in formula (4).
% 
%   kb = KaiserBessel_ImageDomain (x, W, G, beta)
%   
%   Arguments:
%   x - pixel number within image
%   W - Kaiser-Bessel window width in pixel units
%   G - oversampled image size
%   beta - Kaiser-Bessel window shape factor
%   kb - Kaiser-Bessel function in the image domain

% Michal Zarrouk, June 2013.


kb = sinc(sqrt((pi*W.*x./G).^2-beta.^2)/pi);
