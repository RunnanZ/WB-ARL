# Project Configuration Guide

## Project Introduction

This project is implemented based on OpenCV, CUDA, and FFTW3 libraries, mainly including 3D Fourier transform, PSF processing, and back-projection related functions. The main function file of the project is WB-ARL.cpp, serving as the core entry point of the program, it provides three deconvolution implementations to meet various needs, including WB-ARL(CUDA), RLD(CUDA), RLD(C++). Below are detailed environment configuration steps, usage instructions, and function descriptions.

## Core Files & Function Descriptions

1. WB-ARL.cpp：Main function file (program entry point), integrating all core logic for process scheduling and result output.

2. RLD.cu：CUDA parallel optimization core code, containing two key accelerated functions.

3. RLD.cuh：Header file for RLD.cu, declaring CUDA-accelerated related functions for main program calls.

4. WB_back_projector.cpp：Back projector computation module, providing core algorithm support for back projection.

5. WB_back_projector.h：Header file for WB_back_projector.cpp, declaring interface for back projection related functions.

6. FourierTransform.cpp：Non-CUDA-accelerated implementation of Fourier Transform (FFT), providing basic FFT functionality.

7. FourierTransform.h：Header file for FourierTransform.cpp, declaring interfaces for basic FFT related functions.

In WB-ARL.cpp, the project provides three deconvolution implementations with unified function call formats, which can be selected based on requirements:

1. WienerButterworth_ARLDeconvolution: Optimized deconvolution algorithm combining Wiener filtering, Butterworth and CUDA Parallel Acceleration.

2. NativeRLDeconvolution: Traditional deconvolution algorithm combining  CUDA Parallel Acceleration.

3. RLDeconvolution: Traditional deconvolution algorithm for C++.

## Prerequisites

* Operating System: Windows 10/11 (64-bit)

* Development Tool: Visual Studio 2022

* Hardware Requirement: NVIDIA GPU supporting CUDA

## Dependency Download

1. OpenCV

   * Version Requirement: 4.8.0 (other compatible versions can be tried)

   * Download URL: OpenCV Official Website https://opencv.org/releases/

   * Example Installation Path: D:\openCV\opencv

2. CUDA

   * Version Requirement: 12.0 (other compatible versions can be tried)

   * Download URL: NVIDIA CUDA Toolkit https://developer.nvidia.com/cuda-toolkit

   * Example Installation Path: C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0

3. FFTW3

   * Version Requirement: 3.x

   * Download URL: FFTW Official Website https://www.fftw.org/

   * Example Extraction Path: D:\VisualStudio2022\VC\Tools\fftw3

## Project Configuration Steps

1. Create Project

   * Open Visual Studio 2022 and create an "Empty Project"

   * Set project platform to: x64 (Release mode)

   * Add all source files (.cpp) to the project

2. Configure Include Directories

   * Right-click project → Properties → Configuration Properties → VC++ Directories → Include Directories → Edit

   * Add the following paths (replace with actual installation paths):

     * C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0\include

     * D:\openCV\opencv\build\include

     * D:\VisualStudio2022\VC\Tools\fftw3

3. Configure Library Directories

   * Right-click project → Properties → Configuration Properties → VC++ Directories → Library Directories → Edit

   * Add the following paths (replace with actual installation paths):

     * D:\openCV\opencv\build\x64\vc16\lib

     * C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0\lib\x64

     * D:\VisualStudio2022\VC\Tools\fftw3

4. Configure Additional Dependencie

   1. Right-click project → Properties → Configuration Properties → Linker → Input → Additional Dependencies → Edit

   2. Add the following library files：

   3. cufft.lib
      cublas.lib
      cuda.lib
      cudadevrt.lib
      cudart.lib
      cudart_static.lib
      opencv_world480.lib
      OpenCL.lib
      libfftw3-3.lib
      libfftw3f-3.lib
      libfftw3l-3.lib

5. Confirm Configuration

   * Ensure the configuration platform is: x64

   * Ensure the configuration mode is: Release

## Notes

* All paths need to be modified according to the actual installation location

* OpenCV's opencv_world480.lib must correspond to the installed version (e.g., 480 for version 4.8.0)

* FFTW3 requires downloading the precompiled Windows version (64-bit)

* If library linking errors occur, check if the path contains the corresponding .lib files

## IMPORTANT NOTE

Should you have any questions regarding this code and the corresponding results, please contact HeHeng Du (hehengdu@njust.edu.cn).

Reference: 1.Du, Heheng, et al. "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy" 2.Guo, Min, et al. "Rapid image deconvolution and multiview fusion for optical microscopy." Nature biotechnology 38.11 (2020): 1337-1346.
