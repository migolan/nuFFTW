%% nuFFTW command-line inetrface overview
% nuFFTW is an auto-tuning, parallel library for computation of non-uniform
% fast Fourier transforms (nuFFTs).
%
% The command-line interface consists of a several binaries that can be
% called from a terminal:
%
% * |nudft_double|, |nudft_float| - for computation of exact, slow
%   non-uniform discrete Fourier transforms (nuDFTs).
% * |nufft_double|, |nufft_float| - for computation of non-uniform fast
%   Fourier transforms (nuFFTs).
% * |write_nufft_impfile_double|, |write_nufft_impfile_float| - computation
%   of a nuFFT implementations parameters and writing them to a nuFFT
%   implementation file for future use.
% * |nufftw_double|, |nufftw_float| - auto-tune the nuFFT to find the
%   performance-optimal implementation from a space of error-equivalent
%   implementations.
%
% For usage examples, refer to |nuFFT_cli_demo.html| and |nuFFTW_cli_demo.html| .

%% nuDFT (exact (slow) non-uniform discrete Fourier transform)
% * |nudft_double Nthreads direction imagefilename kdatafilename Ndim imagesize trajfilename|
% * |nudft_float  Nthreads direction imagefilename kdatafilename Ndim imagesize trajfilename|

%% nuFFT (non-uniform fast Fourier transform)
% * Specify nuFFT implementation parameters in the command line:
% 
% |nufft_double Nthreads direction imagefilename kdatafilename Ndim imagesize maxerr_power alpha resamplingmethod trajfilename|
%
% |nufft_float  Nthreads direction imagefilename kdatafilename Ndim imagesize maxerr_power alpha resamplingmethod trajfilename|
%
% * Get nuFFT implementation parameters from a pre-written nuFFT
% implementation file:
%
% |nufft_double Nthreads direction imagefilename kdatafilename impfilename|
%
% |nufft_float  Nthreads direction imagefilename kdatafilename impfilename|

%% Write a nuFFT implementation file
% * |write_nufft_impfile_double Nthreads Ndim imagesize maxerr_power alpha resamplingmethod trajfilename impfilename|
% * |write_nufft_impfile_float  Nthreads Ndim imagesize maxerr_power alpha resamplingmethod trajfilename impfilename|

%% nuFFTW - auto-tune the nuFFT
% * Auto-tune and save best implementation to file:
%
% |nufftw_double Nthreads Ndim imagesize maxerr_power alphamin alphamax nalpha resamplingmethod tuning_heuristic trajfilename impfilename|
%
% |nufftw_float  Nthreads Ndim imagesize maxerr_power alphamin alphamax nalpha resamplingmethod tuning_heuristic trajfilename impfilename|
%
% * Auto-tune and just print what the best implementation was:
%
% |nufftw_double Nthreads Ndim imagesize maxerr_power alphamin alphamax nalpha resamplingmethod tuning_heuristic trajfilename|
%
% |nufftw_float  Nthreads Ndim imagesize maxerr_power alphamin alphamax nalpha resamplingmethod tuning_heuristic trajfilename|

%% Argument descriptions
% * |Nthreads| - integer number of threads.
% * |direction| - direction of Fourier transform, either |forward| or |adjoint|.
% * |trajfilename| - name of trajectory data file, which contains consecutive
%   sets of k-space coordinates and dcfs (not sqrt of dcfs!) per each
%   k-space sampling point, written in binary format, in either double or
%   float. (See bottom of this page for an example.)
% * |Ndim| - integer number of dimensions.
% * |imagesize| - set of |Ndim| image sizes in pixels in each dimension
% * |imagefilename| - name of file containing complex uniform data in
%   either double or float precision, as consecutive pairs of real and
%   imaginary values. Data size should correspond to |prod(imagesize)|.
%   (Similar to a .cfl or .cdb file.)
% * |kdatafilename| - name of file containing complex non-uniform data in
%   either double or float precision, as consecutive pairs of real and
%   imaginary values. Data size should correspond to number of non-uniform
%   sampling locations in the file |trajfilename|. (Similar to a .cfl or .cdb file.)
% * |maxerr_power| - maximum aliasing error *(in powers of 10, i.e. for
%   |maxerr = 1e-3|, insert |maxerr_power = -3|)*.
% * |alpha| - grid oversampling ratio.
% * |alphamin| - minimum oversampling ratio for tuning.
% * |alphamax| - maximum oversampling ratio for tuning.
% * |nalpha| - number of oversampling ratios (nuFFT implementations) for tuning.
% * |resamplingmethod| - either |spm| (full precomputation) or |onthefly| (no precomputation).
% * |tuning_heuristic| - either |nufft| (to tune the entire nuFFT
%   computation) or |fftw| (to only tune the FFT stage of the nuFFT).
% * |impfilename| - name of a previously-written nuFFT implementation file.

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