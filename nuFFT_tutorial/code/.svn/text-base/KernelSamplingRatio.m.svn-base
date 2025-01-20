function S = KernelSamplingRatio (maxaliasingamp, alpha, kernel_interpolation_method)
% KernelSamplingRatio - Optimal sampling ratio for the gridding interpolation kernel
% 	as presented by Beatty et al. in formulas (7) and (8).
% 	
%   S = KernelSamplingRatio (maxaliasingamp, alpha, kernel_interpolation_method)
%
% 	Arguments:
% 	maxaliasingamp - maximum aliasing amplitude of gridding, as defined in
%       Beatty et al., formula (3).
% 	alpha - grid oversampling ratio
% 	kernel_interpolation_method - interpolation method of the presampled
%       gridding interpolation kernel: 'NEAREST_NEIGHBOR' or 'LINEAR'

% Michal Zarrouk, June 2013.


switch (kernel_interpolation_method)
    case 'NEAREST_NEIGHBOR'
        S = 0.91/maxaliasingamp/alpha;

    case 'LINEAR'
        S = sqrt(0.37/maxaliasingamp)/alpha;
end
