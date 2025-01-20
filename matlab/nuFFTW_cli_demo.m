%% nuFFTW command-line interface demo
% This demo explains how to use the command-line interface of the nuFFTW library.

%% nuDFT
!../bin/nudft_double 1 adjoint ../data/spiral2d_phantom_imagenudft.cdb ../data/spiral2d_phantom_data.cdb 2 128 128 ../data/spiral2d_phantom_traj.db
imagenudft = reshape(vec2complex(readbin('../data/spiral2d_phantom_imagenudft.cdb','double')),128,128);
%%
% Wow, that's a lot of time for such a small trajectory!

%% nuFFT
!../bin/nufft_double 1 adjoint ../data/spiral2d_phantom_imagenufft.cdb ../data/spiral2d_phantom_data.cdb 2 128 128 -3 2 spm ../data/spiral2d_phantom_traj.db
imagenufft = reshape(vec2complex(readbin('../data/spiral2d_phantom_imagenufft.cdb','double')),128,128);
%%
% That was faster.
%
% Let's check the accuracy:
e = get_image_errors(imagenudft,imagenufft);
%%
% The error between the direct nuFT (nuDFT) and the nuFFT matches the input maximum
% aliasing error with which we created the nuFFT implementation.
%%
figure
subplot(131)
imshowz(imagenudft)
title('adjoint nuDFT')
subplot(132)
imshowz(imagenufft)
title('adjoint nuFFT')
subplot(133)
imshowz(e)
colorbar
title('normalized error')

%% Write nuFFT implementation to file
!../bin/write_nufft_impfile_double 1 2 128 128 -3 2 spm  ../data/spiral2d_phantom_traj.db ../data/spiral2d_phantom.imp
%%
% Now we will perform an adjoint and forward nuFFT using this pre-written
% nuFFT implememntation file
!../bin/nufft_double 1 adjoint ../data/spiral2d_phantom_imagenufft.cdb ../data/spiral2d_phantom_data.cdb ../data/spiral2d_phantom.imp
%%
!../bin/nufft_double 1 forward ../data/spiral2d_phantom_imagenufft.cdb ../data/spiral2d_phantom_data2.cdb ../data/spiral2d_phantom.imp
%%
% Check consistency between the forward and adjoint operators:
trajdata = readbin('../data/spiral2d_phantom_traj.db','double');
Nsamples = length(trajdata)/3;
trajdata = reshape(trajdata,3,Nsamples);
trajectory = trajdata(1:2,:)';
dcf = trajdata(3,:)';

kdata = vec2complex(readbin('../data/spiral2d_phantom_data.cdb','double'));
kdata_dc = kdata .* sqrt(dcf);

imagenufft_vec = vec2complex(readbin('../data/spiral2d_phantom_imagenufft.cdb','double'));

kdata_dc2 = vec2complex(readbin('../data/spiral2d_phantom_data2.cdb','double'));

kdata_dc2'*kdata_dc-imagenufft_vec'*imagenufft_vec

%% nuFFTW - Auto-tune the nuFFT
!../bin/nufftw_double 1 2 128 128 -3 1.2 2 20 spm nufft ../data/spiral2D_phantom_traj.db ../data/spiral2D_phantom_optimal.imp
%%
% Use the pre-written optimal nuFFT implementation file for nuFFT:
!../bin/nufft_double 1 adjoint ../data/spiral2d_phantom_imagenufft.cdb ../data/spiral2D_phantom_data.cdb ../data/spiral2D_phantom_optimal.imp
%%
imagenufft = reshape(vec2complex(readbin('../data/spiral2d_phantom_imagenufft.cdb','double')),128,128);
e = get_image_errors(imagenudft,imagenufft);
figure
subplot(131)
imshowz(imagenudft)
title('adjoint nuDFT')
subplot(132)
imshowz(imagenufft)
title('adjoint nuFFT')
subplot(133)
imshowz(e)
colorbar
title('normalized error')
