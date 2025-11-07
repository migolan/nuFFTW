nuFFTW is an auto-tuning, parallel library for computation of the non-uniform fast Fourier transform (nuFFT).

<p align="center"><img src="nuFFTW ISMRM13 poster.png" alt="ISMRM13 poster" width="600"/></p>

## Installation
For instructions on how to download and install the nuFFTW, please refer to the [installation instructions](INSTALL).

## Usage
The nuFFTW has a command-line interface and a mex interface.
The following demos and overviews of these interfaces are available:

| | nuFFTW mex API | nuFFTW CLI | nuFFT mex API |
|:--------|:--------:|:--------:|:--------:|
| overview    | [✅](matlab/html/nuFFTW_mex_overview.html) | [✅](matlab/html/nuFFTW_cli_overview.html) | |
| html demo   | [✅](matlab/html/nuFFTW_mex_demo.html) | [✅](matlab/html/nuFFTW_cli_demo.html) | [✅](matlab/html/nuFFT_mex_demo.html)|
| matlab demo | [✅](matlab/nuFFTW_mex_demo.m) | [✅](matlab/nuFFTW_cli_demo.m) | [✅](matlab/nuFFT_mex_demo.m) |

Note that the mex interface can be used either in double precision or single but not both at the same time, otherwise MATLAB will likely crash.

## Tutorial
A tutorial on the nuFFT is available [here](nuFFT_tutorial/code/nufft_tutorial).

Please let me know of any problems you encounter during installation or running, so I can track and debug them.
