# WB-ARL: Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution

This repository provides two implementations of the Wiener–Butterworth Accelerated Richardson–Lucy (WB-ARL) deconvolution algorithm for rapid 3D fluorescence microscopy, as described in our paper：for MATLAB Users and  CUDA C++ Users.



## How to use（for matlab users）

The code is tested in MATLAB 2022b (64bit) under the Windows 11 64bit version with an Intel i9-13900KF CPU, NVIDIA GeForce RTX 4060 GPU and 64GB RAM.

1. Unpack the package

2. Include subdirectory in your Matlab path

3. Run the .m files with the prefix of "test" to process the example samples.

a). The raw data from wide-field fluorescence should be placed in the “Raw” folder. Some example data for testing are provided, including 3D stacks of real Chlorella (collected with a 40x/0.75NA objective lens). b). PSF files in the .tif format should be placed in the folder. Readers can also generate simulated ideal PSFs based on wave optics theory and phase space theory.

For details, refer to the WB-ARL(matlab) file.

## How to use（for CUDA C++users）

This project is implemented based on OpenCV, CUDA, and FFTW3 libraries, mainly including 3D Fourier transform, PSF processing, and back-projection related functions. The main function file of the project is WB-ARL.cpp, serving as the core entry point of the program, it provides three deconvolution implementations to meet various needs, including WB-ARL(CUDA), RLD(CUDA), RLD(C++).

1. Operating System: Windows 10/11 (64-bit)

2. Development Tool: Visual Studio 2022

3. Hardware Requirement: NVIDIA GPU supporting CUDA

For details, refer to the WB-ARL(CUDA C++) file.

## Citation

If you use this code please cite the companion paper where the orginal method appeared:

Guo, Min, et al. "Rapid image deconvolution and multiview fusion for optical microscopy." Nature biotechnology 38.11 (2020): 1337-1346. https://doi.org/10.1038/s41587-020-0560-x

## IMPORTANT NOTE

Should you have any questions regarding this code and the corresponding results, please contact HeHeng Du (hehengdu@njust.edu.cn).

Reference: 1.Du, Heheng, et al. "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy" 2.Guo, Min, et al. "Rapid image deconvolution and multiview fusion for optical microscopy." Nature biotechnology 38.11 (2020): 1337-1346.

