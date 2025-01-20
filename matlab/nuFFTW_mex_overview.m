%% nuFFTW mex interface overview
% nuFFTW is an auto-tuning, parallel library for computation of non-uniform
% fast Fourier transforms (nuFFTs).
%
% In order to perform a nuFFT, you should first create a nuFFT
% implementation object.
%
% Once an implementation has been created, the following operations may be
% performed:
%
% * forward transform
% * adjoint transform
% * write implementation to file
% * clear implementation
%
% Other usage options include writing a nuFFT implementation file without
% creating a nuFFT implementation object in MATLAB, and auto-tuning a nuFFT
% to find the performance-optimal implementation from a space of
% error-equivalent implementations.
%
% Note that the mex interface can be used either in double precision or
% single but not both at the same time, otherwise MATLAB will likely crash.
%
% For usage examples, refer to |nuFFT_mex_demo.html| and |nuFFTW_mex_demo.html| .
%
% For complete function reference, type |help nufft_implementation|, |help
% nufft_implementation.(methodname)|, |help write_nufft_impfile|, or |help
% nufftw| (or |doc|).

%% Create a nuFFT implementation object
% * Create implementation from trajectory and square root of density
% compensation factors as MATLAB variables:
%
% |imp = nufft_implementation(datatype, Nthreads, imagesize, maxerr, alpha,
% resamplingmethod, trajectory, sqrtdcf);|
%
% * Create implementation from trajectory and dcfs from a datafile (see
% trajectory datafile format info on the bottom of this page):
%
% |imp = nufft_implementation(datatype, Nthreads, imagesize, maxerr, alpha,
% resamplingmethod, trajfilename);|
%
% * Create implementation from a previously-written nuFFT implementation
% file:
%
% |imp = nufft_implementation(datatype, Nthreads, impfilename);|
 
%% Forward nuFFT
% * |kdata_denscomp = imp.forward(Nthreads, image_data);|
% * |kdata_denscomp = imp*image_data;|

%% Adjoint nuFFT
% Note that the input k-space data should be multiplied by the square root
% of the dcfs.
%
% * |image_data = imp.adjoint(Nthreads, kdata_denscomp);|
% * |image_data = imp'*kdata_denscomp;|

%% Write nuFFT implementation file
% |imp.write_impfile(impfilename);|

%% Clear nuFFT implementation
% |clear imp|

%% Write nuFFT implementation file without creating a nuFFT implementation object in MATLAB
% * Write implementation file from trajectory and square root of density
% compensation factors as MATLAB variables:
%
% |write_nufft_impfile(datatype, Nthreads, imagesize, maxerr, alpha,
% resamplingmethod, trajectory, sqrtdcf, impfilename);|
%
% * Write implementation file from trajectory and dcfs from a datafile (see
% trajectory datafile format info on the bottom of this page):
%
% |write_nufft_impfile(datatype, Nthreads, imagesize, maxerr, alpha,
% resamplingmethod, trajfilename, impfilename);|

%% nuFFTW
% * Auto-tune and write optimal nuFFT implementation to file:
%
% from trajectory and square root of density compensation factors as
% MATLAB variables:
%
% |nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajfilename, impfilename);|
%
% from trajectory and dcfs from a datafile (see trajectory datafile
% format info on the bottom of this page):
%
% |nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajectory, dcf, impfilename);|
%
% * Auto-tune and return optimal nuFFT implementation:
%
% from trajectory and square root of density compensation factors as MATLAB variables:
%
% |imp = nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajfilename);|
%
% from trajectory and dcfs from a datafile (see trajectory datafile format info on the bottom of this page):
%
% |imp = nufftw(datatype, Nthreads, imagesize, maxerr, alphamin, alphamax, nalpha, resamplingmethod, tuning_heuristic, trajectory, dcf);|

%% Argument descriptions
% * |datatype| - either |'double'| or |'float'|, in consistency with the class of the input trajectory data.
% * |Nthreads| - integer number of threads.
% * |imagesize| - array of image size in pixels in each dimension.
% * |maxerr| - maximum aliasing error.
% * |alphamin| - minimum oversampling ratio for tuning.
% * |alphamax| - maximum oversampling ratio for tuning.
% * |nalpha| - number of oversampling ratios (nuFFT implementations) for tuning.
% * |resamplingmethod| - either |'spm'| (full precomputation) or |'onthefly'| (no precomputation).
% * |trajfilename| - name of trajectory data file, which contains consecutive
%   sets of k-space coordinates and dcfs (not sqrt of dcfs!) per each
%   k-space sampling point, written in binary format, in either double or
%   float. (See bottom of this page for an example.)
% * |trajectory| - vector of trajectory coordinates, ordered by sample.
% * |sqrtdcf| - vector of square root of sampling point density
%   compensation factors.
% * |tuning_heuristic| - either |'nufft'| (to tune the entire nuFFT
%   computation) or |'fftw'| (to only tune the FFT stage of the nuFFT).
% * |impfilename| - name of a previously-written nuFFT implementation file.
% * |imp| - nuFFT implementation object.

%% Trajectory datafile format
% A trajectory datafile should contain consecutive sets of k-space
% coordinates and dcfs (not sqrt of dcfs!) per each k-space sampling point,
% written in binary format, in the desired precision (double or float).
%
% For example, a 3-dimensional trajectory's datafile will be formatted as
% following:
%
% |[sample_1_x_coordinate, sample_1_y_coordinate, sample_1_z_coordinate,
% sample_1_dcf, sample_2_x_coordinate, sample_2_y_coordinate,
% sample_2_z_coordinate, sample_2_dcf, ... ]|