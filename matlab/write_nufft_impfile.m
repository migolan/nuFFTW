function write_nufft_impfile(datatype, varargin)
% WRITE_NUFFT_IMPFILE - Compute nuFFT implementation and write to file
%   Call sequences:
% 	  write_nufft_impfile(datatype, Nthreads, imagesize, maxerr, alpha, resamplingmethod, trajectory, sqrtdcf, impfilename);
% 	  write_nufft_impfile(datatype, Nthreads, imagesize, maxerr, alpha, resamplingmethod, trajfilename, impfilename);
%
%   Arguments:
%     datatype - either |'double'| or |'float'|, in consistency
%       with the class of the input trajectory data.
%     Nthreads - integer number of threads.
%     imagesize - array of image size in pixels in each dimension.
%     maxerr - maximum aliasing error.
%     alpha - grid oversampling ratio.
%     resamplingmethod - either |'spm'| (full precomputation) or
%       |'onthefly'| (no precomputation).
%     trajfilename - name of trajectory data file, which
%       contains consecutive sets of k-space coordinates and
%       dcfs (not sqrt of dcfs!) per each k-space sampling
%       point, written in binary format, in either double or float.
%     trajectory - vector of trajectory coordinates, ordered by
%       sample.
%     sqrtdcf - vector of square root of sampling point density
%       compensation factors.
%     impfilename - name of a previously-written nuFFT
%       implementation file.
 
% Michal Zarrouk, July 2013.

if strcmp(datatype,'double')
    nufft_mex_double('write_impfile', varargin{:});
elseif strcmp(datatype,'float')
    nufft_mex_float('write_impfile', varargin{:});
end
