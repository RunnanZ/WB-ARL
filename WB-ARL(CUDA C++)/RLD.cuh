#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <iostream>

#include <cuda_runtime_api.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cufft.h>

//#include <opencv.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
//
#include <complex>
#include "fftw3.h"
//#include "FourierTransform.h"

//#include "Matlab2C.h"
#include <time.h>
#include <chrono>
#include <iostream>
#include <cmath>
#include <string>
#include <math.h>
#include <vector>
#include <ranges>


#ifdef RLD_EXPORTS
#define RLD_API __declspec(dllexport)
#else
#define RLD_API __declspec(dllimport)
#endif


//RLD Ö÷³ÌÐò



extern "C" RLD_API int NativeRLDeconvolution(const std::vector<cv::Mat>&Image_Ori, const std::vector<cv::Mat>&AOTF, std::vector<cv::Mat>&Object_3D_ImageEstimate ,  int NI);
extern "C" RLD_API int WienerButterworth_ARLDeconvolution( std::vector<cv::Mat>&Image_Ori,  std::vector<cv::Mat>&PSF_fp,  std::vector<cv::Mat>&PSF_bp, std::vector<cv::Mat>&Object_3D_ImageEstimate , int NI);//WBRLD


class Matlab2C {

public:

	static void linspace_number(double input_begin, double input_finish, int input_number, cv::Mat& Output_array);  

	static void linspace_interval(double input_begin, double input_finish, double input_interval, cv::Mat& Output_array);  

	static void Repmat(const cv::Mat& input_mat, int repeat_num, std::vector<cv::Mat>& Output_array); 
};

class MatrixOperation {

public:

	static void meshgrid(const cv::Mat input_x, const cv::Mat input_y, cv::Mat& X, cv::Mat& Y);  

	static void exp_ix(const cv::Mat& input_x, cv::Mat& Output_real, cv::Mat& Output_imag);

	static void circshift(cv::Mat& input_array, const cv::Point& move_num);

	static void matSplit(const cv::Mat& input_array, cv::Mat& Output_real, cv::Mat& Output_img); 

	static void matMerge(const cv::Mat& input_real, const cv::Mat& input_imag, cv::Mat& Output_merge);

	static void cMatMultiply(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array); 

	static void Mat_Mul_cMat(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array); 

	static void Mat_Mul_cMat3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output); 

	static void Mat_Mul3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output);


	static void cMatDivide(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array); 

	static void cMat_Div_Mat(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array); 

	static void cMat_Div3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output);


	static void cMat_Div_Mat3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output);


	static void cMatAdd(cv::Mat& input_array1, cv::Mat& input_array2, cv::Mat& Output_array);

	static void Mat_acos(cv::Mat& input_array, cv::Mat& Output_array);

	static void mapminmax(const cv::Mat& input_array, double max_Value, double min_Value, cv::Mat& Output_array);

	static void Nonnegative(cv::Mat& input_array);



};




class TransferFunction {

public:




	static std::vector<cv::Mat> Caculate_F(std::vector<cv::Mat>& rhoxy, std::vector<cv::Mat>& eta, double lambda,
		std::vector<cv::Mat>& rho, double rhos, double rhop);
};






__global__ void kernel_copy(cufftComplex* data1, cufftComplex* data2, int data_len);//real_mul data2 output

__global__ void kernel_copy1(cufftComplex* data1, float* data2, int data_len);//real_mul data2 output

__global__ void kernel_cMat_Div3(cufftComplex* data1, cufftComplex* data2, int data_len);//real div

__global__ void kernel_Mat_Mul3(cufftComplex* data1, cufftComplex* data2, int data_len);//real_mul data2 output

__global__ void kernel_Mat_Mul_cMat3(cufftComplex* data1, cufftComplex* data2, int data_len);

__global__ void kernel_Nonnegative(cufftComplex* data, int data_len);

__global__ void kernel_normalizing(cufftComplex* data, int data_len);//ifft3

__global__ void kernel_normalizing_imag0(cufftComplex* data, int data_len);//ifft3

__global__ void kernel_imag0(cufftComplex* data, int data_len);

__global__ void kernel_fftshift3(cufftComplex* data, cufftComplex* out, int X, int Y, int Z);

__global__ void kernel_ifftshift3(cufftComplex* data, cufftComplex* out, int X, int Y, int Z);