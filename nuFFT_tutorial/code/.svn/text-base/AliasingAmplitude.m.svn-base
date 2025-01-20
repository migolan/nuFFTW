function aliasingamp = AliasingAmplitude(x, W, G, beta)
% AliasingAmplitude - Measure of a pixel's gridding accuracy
%   as defined by Beatty et al. in formula (3).
% 
%   aliasingamp = AliasingAmplitude(x, W, G, beta)
%
%   Arguments:
%   x - pixel number within image
%   W - kernel width
%   G - oversampled image size
%   beta - Kaiser-Bessel kernel shape factor
%   aliasingamp - aliasing amplitude

% Michal Zarrouk, June 2013.


% number of aliasing replications to sum (see formula (3) in Beatty et
% al.). This parameter is actually closely related to maatol from the
% KernelWidth function. In order to get maa = 1e-3, we need around 6
% replications from each side.
sumwidth = 6;

p = 1:sumwidth;
p = [-p p]';
[x1,p] = meshgrid(x,p);

% numerator of formula (3)
argsum = sum((KaiserBessel_ImageDomain(x1+G*p, W, G, beta).^2));

% denominator of formula (3)
a = KaiserBessel_ImageDomain(x, W, G, beta);

aliasingamp = sqrt(argsum)./a;
