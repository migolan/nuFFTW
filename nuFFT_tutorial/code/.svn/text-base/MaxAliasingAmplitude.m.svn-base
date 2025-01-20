function [maxaliasingamp, beta] = MaxAliasingAmplitude(alpha, W, N)
% MaxAliasingAmplitude - Maximum aliasing amplitude of gridding
%   with a given grid oversampling ratio and kernel width, as defined by
%   Beatty et al. in formula (3).
%
%   [maxaliasingamp, beta] = MaxAliasingAmplitude(alpha, W, N)
%
%   Arguments:
%   alpha - grid oversampling ratio
%   W - kernel width in pixel units
%   imagesize - number of pixels in image
%   maxaliasingamp - Maximum aliasing amplitude of gridding
%   beta - optimal Kaiser-Bessel kernel shape factor for gridding with alpha and W

% Michal Zarrouk, June 2013.


% number of pixels in oversampled image
G = OversampledImageSize(alpha,N);

% optimal Kaiser-Bessel kernel shape factor for gridding with alpha and W
beta = KaiserBesselOptimalBeta(W,alpha);

x = -floor(N/2):(ceil(N/2)-1);
aliasingamp = AliasingAmplitude(x, W, G, beta);
maxaliasingamp = max(aliasingamp);
