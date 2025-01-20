function create_nufft_impfile(datatype, varargin)
% CREATE_NUFFT_IMPFILE
% 	create_nufft_impfile(datatype, Nthreads, imagesize, maxerr, alpha, resamplingmethod, trajectory, sqrtdcf, impfilename);
% 	create_nufft_impfile(datatype, Nthreads, imagesize, maxerr, alpha, resamplingmethod, trajfilename, impfilename);

% Michal Zarrouk, July 2013.

if strcmp(datatype,'double')
    nufft_mex_double('create_impfile', varargin{:});
elseif strcmp(datatype,'float')
    nufft_mex_float('create_impfile', varargin{:});
end
