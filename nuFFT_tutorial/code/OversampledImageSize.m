function G = OversampledImageSize(alpha, N)
% OversampledImageSize - Number of pixels in the 1D of the oversampled image
%   
%   G = OversampledImageSize(alpha, N)
%
%   Arguments:
%   alpha - grid oversampling ratio
%   N - number of pixels in each dimension of the image
%   G - number of pixels in each dimension of the oversampled image
%   
%   We enforce an even number so that the phase modulation (which is
%   necessary to avoid fftshifting after the fft) will be just +-1, and not
%   complex (which would require a double amount of storage for the sparse
%   convolution matrix, and complex multiplication).

% Michal Zarrouk, June 2013.


G = ceil(alpha*N);
if mod(G,2) == 1
    G = G+1;
end
