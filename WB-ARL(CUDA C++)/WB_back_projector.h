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
#include <omp.h>
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
//class PSFOperation {
//public:
//	// ... 其他函数声明 ...
//
//	// PSF翻转函数
//	static void flipPSF(const std::vector<cv::Mat>& inPSF, std::vector<cv::Mat>& outPSF);
//};
//
void Generate_back_projector(
	std::vector<cv::Mat>& PSF,      // forward projector PSF
	const std::string& bp_type,               // back projector
	double alpha,                             
	double beta,                              
	int n_order,                                   
	const std::vector<double>& iRes,          
	
	std::vector<cv::Mat>& PSF_bp);             




void flipPSF(std::vector<cv::Mat>& inPSF, std::vector<cv::Mat>& outPSF);
void Abs_complex(const std::vector<cv::Mat>& inputMats, std::vector<cv::Mat>& outputMats);
void Abs_complex(const cv::Mat& inputMat, cv::Mat& outputMat);//overload
double Max(const std::vector<cv::Mat>& mats);
cv::Mat rowMax(const cv::Mat& input);

std::vector<cv::Mat> calculatePSF_bp(
	 std::vector<cv::Mat>& flippedPSF,
	const std::string& bp_type,
	double alpha,
	double beta,
	int n_order,
	const std::vector<double>& iRes);

