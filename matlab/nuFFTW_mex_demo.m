%% nuFFTW mex interface demo
% This example demonstrates how to use the mex interface of the nuFFTW.
%
% Before reading this, you should familiarize yourself with the nuFFTW mex
% interface by reading the nuFFTW mex interface overview
% (nuFFTW_mex_overview.html) and the nuFFT mex interface demo
% (nuFFT_mex_demo.html).
%
% nuFFTW is an auto-tuning, parallel library for computation of non-uniform
% fast Fourier transforms (nuFFTs). It provides the possibility to create
% an nuFFT implementation with the highest performance for a given level of
% gridding accuracy. To that end, we define a space of error-equivalent
% nuFFT implementations with different oversampling ratios (alpha).
% This space is specified by a minimum and maximum oversampling ratio
% (\alpha), and the number of implementations to be checked.

addpath ../bin
load ../data/spiral2d_phantom.mat
trajectory = kcoord';
trajectory = trajectory(:);
sqrtdcf = sqrt(dcf);
Nthreads = 1;
imagesize = [128 128];
maxaliasingerror = 1e-3;
alphamin = 1.2;
alphamax = 2;
Nalpha = 20;
resampling_method = 'spm';
tuning_heuristic = 'nufft';

%%
% Tune the nuFFT implementation space:
nuFFT_imp = nufftw('double', Nthreads, imagesize, maxaliasingerror, alphamin, alphamax, Nalpha, resampling_method, tuning_heuristic, trajectory, sqrtdcf);

%%
% Note that the implementation counting is zero-based.
%
% This implementation can now be used for nuFFTs:
kdata_dc = kdata.*sqrtdcf;
imagenufft = nuFFT_imp.adjoint(1,kdata_dc);
e = get_image_errors(imagenudft,imagenufft);
figure
imshowz(imagenufft)
title('adjoint nuFFT')

