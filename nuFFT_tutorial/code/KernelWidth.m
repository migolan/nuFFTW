function [W, beta] = KernelWidth(maa, alpha, N)
% KernelWidth - Width of gridding interpolation kernel
%   that corresponds to a given grid oversampling ratio and maximum
%   aliasing amplitude, as defined by Beatty et al. in formula (3) (see
%   fig. 3).
% 
%   [W, beta] = KernelWidth(maa, alpha, N)
%
%   Arguments:
%   maa - maximum aliasing amplitude of gridding
%   alpha - grid oversampling ratio
%   N - number of pixels in each dimension of the image
%   W - width of gridding interpolation kernel, in pixel units
%   beta - optimal Kaiser-Bessel kernel shape factor that corresponds to alpha and W

% Michal Zarrouk, June 2013.


% define the search range for the kernel width
Wmin = 1;
Wmax = 10;

% percentage of error tolerance allowed for deviation from maa
maatol = 1e-3;

% maximum aliasing amplitude of current estimate of kernel width
curmaa = realmax;

% Binary search for the kernel width that corresponds to alpha and mae
while abs((curmaa - maa)/maa) > maatol
    W = (Wmax + Wmin)/2;
    [curmaa, beta] = MaxAliasingAmplitude(alpha, W, N);
    if curmaa > maa
        Wmin = W;
    else
        Wmax = W;
    end
end
