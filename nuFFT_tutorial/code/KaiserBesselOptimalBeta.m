function beta = KaiserBesselOptimalBeta (W, alpha)
% KaiserBesselOptimalBeta - Optimal shape factor for gridding
%   with the Kaiser-Bessel kernel, as defined by Beatty et al. in formula (5).
% 
%   beta = KaiserBesselOptimalBeta (W, alpha)
%   
%   Arguments:
%   W - width of the gridding interpolation kernel
%   alpha - grid oversampling ratio
%   beta - Kaiser-Bessel kernel shape factor

% Michal Zarrouk, June 2013.


beta = pi * sqrt((W./alpha.*(alpha-0.5)).^2-0.8);
