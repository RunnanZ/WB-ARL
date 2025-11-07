#include "WB_back_projector.h"
#include <opencv2/opencv.hpp>
#include <iostream>  
#include "RLD.cuh"
#include "FourierTransform.h"
#include <omp.h>
using namespace std;
void Generate_back_projector(
	std::vector<cv::Mat>& PSF,     
	const std::string& bp_type,              
	double alpha,                             
	double beta,                              
	int n_order,                                    
	const std::vector<double>& iRes,          
	std::vector<cv::Mat>& PSF_bp) {



void flipPSF(std::vector<cv::Mat>& inPSF, std::vector<cv::Mat>& outPSF) {

	const int Sz = inPSF.size();          // Sz 
	const int Sx = inPSF[0].rows;         // Sx 
	const int Sy = inPSF[0].cols;         // Sy 


	outPSF.resize(Sz);
	for (int k = 0; k < Sz; k++) {
		outPSF[k] = cv::Mat::zeros(Sx, Sy, inPSF[0].type());
	}

	for (int i = 0; i < Sx; i++) {       
		 int flipped_i = Sx - 1 - i;  

		for (int j = 0; j < Sy; j++) {    
			 int flipped_j = Sy - 1 - j; 

			for (int k = 0; k < Sz; k++) { 
				 int flipped_k = Sz - 1 - k;



					double val = inPSF[flipped_k].at<double>(flipped_i, flipped_j);
					//cout<< val <<endl;
					outPSF[k].at<double>(i, j) = val;

					//// 
					//if (i == 379 && j == 380 && k == 38) {
					//	std::cout << "Debug [0,0,0]: input value at (" << flipped_i << "," << flipped_j << "," << flipped_k << ") = " << val << std::endl;
					//}
					//if (i == Sx - 1 && j == Sy - 1 && k == Sz - 1) {
					//	std::cout << "Debug [end]: input value at (" << flipped_i << "," << flipped_j << "," << flipped_k << ") = " << val << std::endl;
					//}

			}
		}
	}
}


void Abs_complex(const std::vector<cv::Mat>& inputMats,  std::vector<cv::Mat>& outputMats) {
	
	outputMats.clear();
	outputMats.reserve(inputMats.size());

	for (const auto& inputMat : inputMats) {
		
		CV_Assert(inputMat.type() == CV_64FC2);

		
		std::vector<cv::Mat> channels(2);
		cv::split(inputMat, channels);
		
		cv::Mat reSq, imSq, magnitudeSq;
		cv::multiply(channels[0], channels[0], reSq);
		cv::multiply(channels[1], channels[1], imSq);
		cv::add(reSq, imSq, magnitudeSq);

		cv::Mat magnitude;
		cv::sqrt(magnitudeSq, magnitude); 

		outputMats.emplace_back();
		magnitude.convertTo(outputMats.back(), CV_64FC1);

	}
}

void Abs_complex(const cv::Mat& inputMat, cv::Mat& outputMat) { //
	
	if (inputMat.channels() == 2) {
		CV_Assert(inputMat.type() == CV_64FC2);


		std::vector<cv::Mat> channels(2);
		cv::split(inputMat, channels);
		cv::Mat realPart = channels[0];
		cv::Mat imagPart = channels[1];


		cv::Mat reSq, imSq, magnitudeSq;
		cv::multiply(realPart, realPart, reSq);
		cv::multiply(imagPart, imagPart, imSq);
		cv::add(reSq, imSq, magnitudeSq);


		cv::sqrt(magnitudeSq, outputMat);

		outputMat.convertTo(outputMat, CV_64FC1);
	}
	
}


double Max(const std::vector<cv::Mat>& mats) {
	double globalMax = -DBL_MAX; 

	for (const auto& mat : mats) {

		CV_Assert(mat.channels() == 1);
		double minVal, maxVal;
		cv::minMaxLoc(mat, &minVal, &maxVal);
		if (maxVal > globalMax) {
			globalMax = maxVal;
		}
	}

	return globalMax;
}

cv::Mat rowMax(const cv::Mat& input) {

	cv::Mat maxValues;
	cv::reduce(input, maxValues, 1, cv::REDUCE_MAX);
	return maxValues;
}





std::vector<cv::Mat> calculatePSF_bp(
	std::vector<cv::Mat>& flippedPSF,
	const std::string& bp_type ,
	double alpha ,
	double beta ,
	int n_order ,
	const std::vector<double>& iRes )
{

	std::cout << "caculate PSF_bp" << std::endl;


	const int size_z = static_cast<int>(flippedPSF.size());
	const int size_y = flippedPSF[0].rows;  
	const int size_x = flippedPSF[0].cols;  

	const double Scx = (static_cast<double>(size_x) + 1.0) / 2.0;
	const double Scy = (static_cast<double>(size_y) + 1.0) / 2.0;
	const double Scz = (static_cast<double>(size_z) + 1.0) / 2.0;

	int Sox = static_cast<int>(std::round((static_cast<double>(size_x) + 1.0) / 2.0));
	int Soy = static_cast<int>(std::round((static_cast<double>(size_y) + 1.0) / 2.0));
	int Soz = static_cast<int>(std::round((static_cast<double>(size_z) + 1.0) / 2.0));
	cout << "Sox=" << Sox << endl;
	cout << "Soy=" << Soy << endl;
	cout << "Soz=" << Soz << endl;
	int Nz = size_z;


	std::vector<cv::Mat> PSF_bp(Nz);
	std::vector<cv::Mat> OTF_bp(Nz);

	// Generate_back_projector
	std::vector<cv::Mat> OTF_flip(Nz), OTF_abs(Nz);


	FourierTransform::ifftshift3(flippedPSF);//		
	FourierTransform::fft3(flippedPSF, OTF_flip, 2);//OTF_flip = fftn(ifftshift(flippedPSF));  
	//cv::Mat OTF_flip37 = OTF_flip[37];// OTF_flip test对

	Abs_complex(OTF_flip, OTF_abs);
	//cv::Mat OTF_abs37 = OTF_abs[37];// OTF_abs test对
	//cv::Mat OTF_abs0 = OTF_abs[0];
	FourierTransform::fftshift3(OTF_abs);//OTF_abs = fftshift(abs(OTF_flip));

	double OTFmax = Max(OTF_abs);//对 OTFmax==M  matlab：OTFmax = max(OTF_abs(:));
	cout << OTFmax << endl;





	const double px = 1.0 / static_cast<double>(size_x);
	const double py = 1.0 / static_cast<double>(size_y);
	const double pz = 1.0 / static_cast<double>(size_z);

	const double tx = 1.0 / (iRes[0] * px);
	const double ty = 1.0 / (iRes[1] * py);
	const double tz = 1.0 / (iRes[2] * pz);


	std::cout << "px: " << px << " py: " << py << " pz: " << pz << std::endl;
	std::cout << "tx: " << tx << " ty: " << ty << " tz: " << tz << std::endl;

	if (bp_type == "wiener") {
		// OTF_flip_norm
		std::vector<cv::Mat> OTF_flip_norm(Nz);
		for (int zz = 0; zz < Nz; zz++) {

			OTF_flip[zz].convertTo(OTF_flip_norm[zz], CV_64FC2);
			OTF_flip_norm[zz] = OTF_flip_norm[zz] / OTFmax;
		}


		for (int zz = 0; zz < Nz; zz++) {

			cv::Mat currentOTF = OTF_flip_norm[zz];

			std::vector<cv::Mat> channels;
			cv::split(currentOTF, channels);
			cv::Mat realPart = channels[0];
			cv::Mat imagPart = channels[1];


			cv::Mat absSquared;
			cv::magnitude(realPart, imagPart, absSquared); 
			cv::pow(absSquared, 2, absSquared);            

			// |OTF|² + alpha
			cv::Mat denominator = absSquared + alpha;
			denominator.convertTo(denominator, currentOTF.type());
			//// realPart / denominator
			//cv::Mat realResult;
			//cv::divide(realPart, denominator, realResult);

			//// imagPart / denominator
			//cv::Mat imagResult;
			//cv::divide(imagPart, denominator, imagResult);

			////
			//std::vector<cv::Mat> resultChannels = { realResult, imagResult };
			//cv::merge(resultChannels, OTF_bp[zz]);
			//////cv::Mat complexResult;
			//////cv::divide(currentOTF, denominator, complexResult); // 直接使用OpenCV复数除法


			//complexResult.copyTo(OTF_bp[zz]);
			// 
			cv::Mat denominator2C;
			std::vector<cv::Mat> denomChannels = { denominator, denominator };
			cv::merge(denomChannels, denominator2C); 

			currentOTF.convertTo(currentOTF, CV_64FC2);
			denominator2C.convertTo(denominator2C, CV_64FC2);

			cv::Mat complexResult;
			cv::divide(currentOTF, denominator2C, complexResult);

			complexResult.copyTo(OTF_bp[zz]);
		}

		FourierTransform::ifft3(OTF_bp, PSF_bp, 1);																//  PSF_bp = fftshift(real(ifftn(OTF_bp)));
		FourierTransform::fftshift3(PSF_bp);

	}
	else if (bp_type == "wiener-butterworth") {

		//OTF_flip_norm 
		std::vector<cv::Mat> OTF_flip_norm(Nz);
		for (int zz = 0; zz < Nz; zz++) {
			OTF_flip[zz].convertTo(OTF_flip_norm[zz], CV_64FC2);
			OTF_flip_norm[zz] = OTF_flip_norm[zz] / OTFmax;
		}


		std::vector<cv::Mat> OTF_wiener(Nz);

		// OTF_wiener 
		for (int zz = 0; zz < Nz; zz++) {

			cv::Mat currentOTF = OTF_flip_norm[zz];


			std::vector<cv::Mat> channels;
			cv::split(currentOTF, channels);
			cv::Mat realPart = channels[0];
			cv::Mat imagPart = channels[1];

			// |OTF|²（模的平方）
			cv::Mat absSquared;
			cv::magnitude(realPart, imagPart, absSquared); 
			cv::pow(absSquared, 2, absSquared);          

			//  |OTF|² + alpha
			cv::Mat denominator = absSquared + alpha;
			denominator.convertTo(denominator, currentOTF.type());

			cv::Mat denominator2C;
			std::vector<cv::Mat> denomChannels = { denominator, denominator };
			cv::merge(denomChannels, denominator2C); // 


			currentOTF.convertTo(currentOTF, CV_64FC2);
			denominator2C.convertTo(denominator2C, CV_64FC2);


			cv::Mat complexResult;
			cv::divide(currentOTF, denominator2C, complexResult);

	
			complexResult.copyTo(OTF_wiener[zz]);
		}

		std::vector<cv::Mat> OTF_Wiener_abs(Nz);

		Abs_complex(OTF_wiener, OTF_Wiener_abs);

		FourierTransform::fftshift3(OTF_Wiener_abs);//OTF_abs = fftshift(abs(OTF_flip));		
		//cv::Mat OTF_Wiener_abs0 = OTF_Wiener_abs[0];//OTF_Wiener_abs  test对
		//cv::Mat OTF_Wiener_abs37 = OTF_Wiener_abs[37];//	

		cv::Mat tplane = OTF_Wiener_abs[Soz - 1];//
		//tplane=cv::abs(tplane); 
		//Abs_complex(tplane, tplane);
		cv::Mat tline;
		tline = rowMax(tplane);

		// to1
		int to1 = std::max(std::round(Scx - tx), 1.0); 

		//  to2
		int to2 = std::min(std::round(Scx + tx), static_cast<double>(size_x));   
		std::cout << "to1: " << to1 << " to2: " << to2 << std::endl;

		double value1 = tline.at<double>(to1 - 1);
		double value2 = tline.at<double>(to2 - 1);
		double beta_wienerx = (value1 + value2) / 2.0;
		std::cout << "value1: " << value1 << " value2: " << value2 << " beta_wienerx: " << beta_wienerx << std::endl;
		double ee = beta_wienerx / (beta * beta) - 1.0;
		std::cout << "ee: " << ee << std::endl;

		std::vector<cv::Mat> mask(Nz);
		for (int k = 0; k < Nz; ++k) {
			mask[k] = cv::Mat::zeros(size_y, size_x, CV_64FC1); 
		}

		cv::Mat tline1111 = tline;
	
		for (int i = 0; i < size_x; ++i) {
			for (int j = 0; j < size_y; ++j) {
				for (int k = 0; k < size_z; ++k) {
					
					double dx = (i + 1 - Scx) / tx;
					double dy = (j + 1 - Scy) / ty;
					double dz = (k + 1 - Scz) / tz;
					double w = dx * dx + dy * dy + dz * dz;

					double mask_value = 1.0 / sqrt(1.0 + ee * pow(w, n_order));

					mask[k].at<double>(j, i) = mask_value;
				}
			}
		}
		//cv::Mat mask10 = mask[10];//mask  test
		//cv::Mat mask0 = mask[0];
		//cv::Mat mask20 = mask[20];
		FourierTransform::ifftshift3(mask);
		//  matlabcode			OTF_bp = mask.*OTF_Wiener; % final OTF_bp cutfoff gain : beta
		MatrixOperation::Mat_Mul_cMat3(mask, OTF_wiener, OTF_bp);
		FourierTransform::ifft3(OTF_bp, PSF_bp, 1);																//  PSF_bp = fftshift(real(ifftn(OTF_bp)));
		FourierTransform::fftshift3(PSF_bp);//PSF_bp  test
		std::cout << "complete caculate PSF_bp" << std::endl;
	}

	return PSF_bp;
}