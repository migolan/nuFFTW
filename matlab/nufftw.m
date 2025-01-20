function varargout = nufftw(datatype, varargin)
% nuFFTW - Auto-tune the nuFFT
%   Call sequence:
%     Auto-tune and write optimal nuFFT implementation to file:
%       nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajfilename, impfilename);
%       nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajectory, sqrtdcf, impfilename);
%     Auto-tune and return optimal nuFFT implementation:
%       imp = nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajfilename);
%       imp = nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajectory, sqrtdcf);
%   
%   Arguments:
%     datatype - either |'double'| or |'float'|, in consistency with the class
%       of the input trajectory data.
%     Nthreads - integer number of threads.
%     imagesize - array of image size in pixels in each dimension.
%     maxerr - maximum aliasing error.
%     alphamin - minimum oversampling ratio for tuning.
%     alphamax - maximum oversampling ratio for tuning.
%     nalpha - number of oversampling ratios (nuFFT implementations) for tuning.
%     resamplingmethod - either |'spm'| (full precomputation) or |'onthefly'|
%       (no precomputation).
%     trajfilename - name of trajectory data file, which contains
%       consecutive sets of k-space coordinates and dcfs (not sqrt of
%       dcfs!) per each k-space sampling point, written in binary format,
%       in either double or float.
%     trajectory - vector of trajectory coordinates, ordered by sample.
%     sqrtdcf - vector of square root of sampling point density
%       compensation factors.
%     tuning_heuristic - either |'nufft'| (to tune the entire nuFFT
%       computation) or |'fftw'| (to only tune the FFT stage of the nuFFT).
%     impfilename - name of a previously-written nuFFT implementation file.
%     imp - nuFFT implementation object.

% Michal Zarrouk, July 2013.

if strcmp(datatype,'double')
    if nargout == 0
        nufftw_mex_double(varargin{:});
        varargout = {};
    else
        [varargout{1:nargout}] = nufftw_mex_double(varargin{:});
    end
elseif strcmp(this.datatype,'float')
    if nargout == 0
        nufftw_mex_float(varargin{:});
        varargout = {};
    else
        varargout{:} = nufftw_mex_float(varargin{:});
    end
end

if nargout > 0
    varargout{1} = nufft_implementation(datatype, varargout{1});
end

