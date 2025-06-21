# WB-ARL
WB-ARL
Title:	 Matlab code for "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy"
Version: 1.0 
Copyright: 2025, HEHENG DU, RUNNAN ZHANG, AND CHAO ZUO

Edited based on the reference [1][2].

Matlab code for "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy"
==========================================================

This package contains the implementations of the algorithms described in the paper, including Wiener-Butterworth back-projection generation and reconstruction algorithms for 3D fluorescence deconvolution.  which are described in detail in the paper:  HEHENG DU, RUNNAN ZHANG, and CHAO ZUO etc, "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy".

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) in an academic publication.

For algorithmic details, please refer to our paper.

----------------
How to use
----------------
The code is tested in MATLAB 2022b (64bit) under the  Windows 11 64bit version with an Intel i9-13900KF CPU, NVIDIA GeForce RTX 4060 GPU and 64GB RAM.

1. Unpack the package
2. Include subdirectory in your Matlab path
3. Run the .m files with the prefix of "test" to process the example samples.

a). The raw data from wide-field fluorescence should be placed in the “Raw” folder. Some example data for testing are provided, including 3D stacks of real Chlorella (collected with a 40x/0.75NA objective lens).
b). PSF files in the .tif format should be placed in the folder. Readers can also generate simulated ideal PSFs based on wave optics theory and phase space theory.



----------------
Main modules description
----------------
1. WB_back_projector.m: Generates Wiener-Butterworth (WB) back-projector based on forward projector (PSF).
2. test.m: Performs 3D reconstruction of WB back-projection deconvolution on a 3D fluorescence stack.


----------------
Citation 
---------------- 
If you use this code please cite the companion paper where the orginal method appeared:

Guo, Min, et al. "Rapid image deconvolution and multiview fusion for optical microscopy." Nature biotechnology 38.11 (2020): 1337-1346. https://doi.org/10.1038/s41587-020-0560-x

----------------
IMPORTANT NOTE 
---------------- 
Should you have any questions regarding this code and the corresponding results, please contact Heheng Du (hehengdu@njust.edu.cn), Runnan Zhang (runnanzhang@njust.edu.cn), or  Chao Zuo (zuochao@njust.edu.cn).



Reference:
1.Du, Heheng, et al. "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy"

2.Guo, Min, et al. "Rapid image deconvolution and multiview fusion for optical microscopy." Nature biotechnology 38.11 (2020): 1337-1346. 

