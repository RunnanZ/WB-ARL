#pragma once
#include <iostream>
#include <opencv2/opencv.hpp>
#include "RLD.cuh"
#include "fftw3.h"

class FourierTransform {

public:

	static void fft2(cv::Mat& input_array, cv::Mat& Output_array);  //fft2

	static void fft3(std::vector<cv::Mat>& input_img, std::vector<cv::Mat>& output, int outputn = 2);  //fft3

	static void ifft2(cv::Mat& input_array, cv::Mat& Output_array);//ifft

	static void ifft3(std::vector<cv::Mat>& input_img, std::vector<cv::Mat>& output, int outputn = 2);  //ifft3

	static void fftshift(cv::Mat& input_array);  // £»

	static void fftshift3(std::vector<cv::Mat>& input_img);

	static void ifftshift(cv::Mat& input_array);  // £»

	static void ifftshift3(std::vector<cv::Mat>& input_img);


};
