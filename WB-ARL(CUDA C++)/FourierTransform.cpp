#include "FourierTransform.h"

void FourierTransform::fft2(cv::Mat& input_array, cv::Mat& Output_array)  //fft2                                                          
{                                                                         
	int m = cv::getOptimalDFTSize(input_array.rows);                         
	int n = cv::getOptimalDFTSize(input_array.cols);

	cv::Mat expand_array, dft_array;
	cv::copyMakeBorder(input_array, expand_array, 0, m - input_array.rows, 0, n - input_array.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));

	if (input_array.channels() == 1) {  
		cv::Mat planes[] = { cv::Mat_<double>(expand_array), cv::Mat::zeros(expand_array.size(), CV_64FC1) };
		cv::Mat complexI;
		merge(planes, 2, complexI);
		dft(complexI, dft_array, cv::DFT_COMPLEX_OUTPUT);
	}
	else if (input_array.channels() == 2) { 
		dft(expand_array, dft_array, cv::DFT_COMPLEX_OUTPUT);
	}

	if (m != 0 || n != 0) {  
		cv::resize(dft_array, Output_array, cv::Size(m, n), 0, 0, cv::INTER_LINEAR); 
	}
	else Output_array = dft_array.clone();
}

void FourierTransform::ifft2(cv::Mat& input_array, cv::Mat& Output_array)//ifft
{
	int m = cv::getOptimalDFTSize(input_array.rows);
	int n = cv::getOptimalDFTSize(input_array.cols);

	cv::Mat expand_array, idft_array;
	cv::copyMakeBorder(input_array, expand_array, 0, m - input_array.rows, 0, n - input_array.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));

	if (input_array.channels() == 1) {  
		cv::Mat planes[] = { cv::Mat_<double>(expand_array), cv::Mat::zeros(expand_array.size(), CV_64FC1) };
		cv::Mat complexI;
		merge(planes, 2, complexI);
		idft(complexI, idft_array, cv::DFT_SCALE + cv::DFT_COMPLEX_OUTPUT);
	}
	else if (input_array.channels() == 2) { 
		idft(expand_array, idft_array, cv::DFT_SCALE + cv::DFT_COMPLEX_OUTPUT);
	}

	if (m != 0 || n != 0) {  
		cv::resize(idft_array, Output_array, cv::Size(m, n), 0, 0, cv::INTER_LINEAR); 
	}
	else Output_array = idft_array.clone();
}

void FourierTransform::fftshift(cv::Mat& input_array)  
{
	cv::Size sz = input_array.size();
	cv::Point move_num(0, 0);
	move_num.x = (int)floor(sz.width / 2.0); // floor
	move_num.y = (int)floor(sz.height / 2.0);
	MatrixOperation::circshift(input_array, move_num);
}

void FourierTransform::ifftshift(cv::Mat& input_array)  
{
	cv::Size sz = input_array.size();
	cv::Point move_num(0, 0);
	move_num.x = -(int)floor(sz.width / 2.0);
	move_num.y = -(int)floor(sz.height / 2.0);
	MatrixOperation::circshift(input_array, move_num);
}



void FourierTransform::fftshift3(std::vector<cv::Mat>& input_img)
{

	int Nz = input_img.size();
	int Nx = input_img[0].rows;
	int Ny = input_img[0].cols;

	std::vector<cv::Mat>  input_img_temp(Nz);

	for (int zz = 0; zz < Nz; zz++) {
		input_img_temp[zz] = input_img[zz].clone();
		FourierTransform::fftshift(input_img_temp[zz]);
	}

	cv::Mat Zshift = cv::Mat::zeros(cv::Size(1, Nz), CV_64FC1);
	for (int zz = 0; zz < Nz; zz++)
	{
		Zshift.at<double>(zz) = zz;
		//std::cout << num << std::endl;
	}
	FourierTransform::fftshift(Zshift);

	for (int zz = 0; zz < Nz; zz++) {

		input_img[zz] = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);
		int Zshif = Zshift.at<double>(zz);
		input_img[zz] = input_img_temp[Zshif];

	}
};

void FourierTransform::ifftshift3(std::vector<cv::Mat>& input_img)
{

	int Nz = input_img.size();
	int Nx = input_img[0].rows;
	int Ny = input_img[0].cols;

	std::vector<cv::Mat>  input_img_temp(Nz);
	for (int zz = 0; zz < Nz; zz++) {
		input_img_temp[zz] = input_img[zz].clone();
		FourierTransform::ifftshift(input_img_temp[zz]);
	}
	cv::Mat Zshift = cv::Mat::zeros(cv::Size(1, Nz), CV_64FC1);
	for (int zz = 0; zz < Nz; zz++)
	{
		Zshift.at<double>(zz) = zz;
		//std::cout << num << std::endl;
	}
	FourierTransform::ifftshift(Zshift);

	for (int zz = 0; zz < Nz; zz++) {

		input_img[zz] = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);
		int Zshif = Zshift.at<double>(zz);
		input_img[zz] = input_img_temp[Zshif];

	}
};


void FourierTransform::fft3(std::vector<cv::Mat>& input_img, std::vector<cv::Mat>& output, int outputn)
{
	int Nz = input_img.size();
	int Nx = input_img[0].cols;
	int Ny = input_img[0].rows;
	int channel = input_img[0].channels();

	fftw_complex* input_array;
	fftw_complex* output_array;

	input_array = (fftw_complex*)fftw_malloc(Nz * Ny * Nx * sizeof(fftw_complex));//
	output_array = (fftw_complex*)fftw_malloc(Nz * Ny * Nx * sizeof(fftw_complex));

	if (channel == 1) {
		for (int zz = 0; zz < Nz; zz++) {

			for (int yy = 0; yy < Ny; yy++)
			{
				for (int xx = 0; xx < Nx; xx++)
				{
					double pv = input_img[zz].at<double>(yy, xx);
					input_array[(zz * Ny + yy) * Nx + xx][0] = pv;//
					input_array[(zz * Ny + yy) * Nx + xx][1] = 0;
				}
			}

		}
	}

	else if (channel == 2)
	{
		cv::Mat real, imag;

		for (int zz = 0; zz < Nz; zz++) {

			MatrixOperation::matSplit(input_img[zz], real, imag);

			for (int yy = 0; yy < Ny; yy++)
			{
				for (int xx = 0; xx < Nx; xx++)
				{
					double pv_real = real.at<double>(yy, xx);
					double pv_imag = imag.at<double>(yy, xx);

					input_array[(zz * Ny + yy) * Nx + xx][0] = pv_real;
					input_array[(zz * Ny + yy) * Nx + xx][1] = pv_imag;
				}
			}

		}


	}

	fftw_plan forward = fftw_plan_dft_3d(Nz, Ny, Nx, input_array, output_array, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(forward);
	fftw_destroy_plan(forward);
	fftw_cleanup();
	if (outputn == 2)
	{
		cv::Mat output_real, output_imag;

		for (int zz = 0; zz < Nz; ++zz)
		{
			output_real = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);
			output_imag = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);

			//std::cout << zz << std::endl;
			for (int yy = 0; yy < Ny; ++yy)
			{
				for (int xx = 0; xx < Nx; ++xx)
				{
					output_real.at<double>(yy, xx) = output_array[(zz * Ny + yy) * Nx + xx][0];
					output_imag.at<double>(yy, xx) = output_array[(zz * Ny + yy) * Nx + xx][1];
				}
			}
			MatrixOperation::matMerge(output_real, output_imag, output[zz]);
		}
	}

	else if (outputn == 1)
	{

		for (int zz = 0; zz < Nz; ++zz)
		{
			cv::Mat output_real = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);

			//std::cout << zz << std::endl;
			for (int yy = 0; yy < Ny; ++yy)
			{
				for (int xx = 0; xx < Nx; ++xx)
				{
					output_real.at<double>(yy, xx) = output_array[(zz * Ny + yy) * Nx + xx][0];
				}
			}

			output[zz] = output_real;
		}

	}

	fftw_free(input_array);
	fftw_free(output_array);
}
//

void FourierTransform::ifft3(std::vector<cv::Mat>& input_img, std::vector<cv::Mat>& output, int outputn)
{
	int Nz = input_img.size();
	int Nx = input_img[0].cols;
	int Ny = input_img[0].rows;
	int channel = input_img[0].channels();
	int NN = Nz * Nx * Ny;

	fftw_complex* input_array;
	fftw_complex* output_array;

	input_array = (fftw_complex*)fftw_malloc(Nz * Ny * Nx * sizeof(fftw_complex));
	output_array = (fftw_complex*)fftw_malloc(Nz * Ny * Nx * sizeof(fftw_complex));

	if (channel == 1) {
		for (int zz = 0; zz < Nz; zz++) {

			for (int yy = 0; yy < Ny; yy++)
			{
				for (int xx = 0; xx < Nx; xx++)
				{
					double pv = input_img[zz].at<double>(yy, xx);
					input_array[(zz * Ny + yy) * Nx + xx][0] = pv;
					input_array[(zz * Ny + yy) * Nx + xx][1] = 0;
				}
			}

		}
	}

	else if (channel == 2)
	{
		cv::Mat real, imag;

		for (int zz = 0; zz < Nz; zz++) {

			MatrixOperation::matSplit(input_img[zz], real, imag);

			for (int yy = 0; yy < Ny; yy++)
			{
				for (int xx = 0; xx < Nx; xx++)
				{
					double pv_real = real.at<double>(yy, xx);
					double pv_imag = imag.at<double>(yy, xx);

					input_array[(zz * Ny + yy) * Nx + xx][0] = pv_real;
					input_array[(zz * Ny + yy) * Nx + xx][1] = pv_imag;
				}
			}
		}
	}


	fftw_plan forward = fftw_plan_dft_3d(Nz, Ny, Nx, input_array, output_array, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(forward);
	fftw_destroy_plan(forward);
	fftw_cleanup();

	if (outputn == 2)
	{
		cv::Mat output_real, output_imag;

		for (int zz = 0; zz < Nz; ++zz)
		{
			output_real = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);
			output_imag = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);

			//std::cout << zz << std::endl;
			for (int yy = 0; yy < Ny; ++yy)
			{
				for (int xx = 0; xx < Nx; ++xx)
				{
					output_real.at<double>(yy, xx) = output_array[(zz * Ny + yy) * Nx + xx][0] / (Nx * Ny * Nz);
					output_imag.at<double>(yy, xx) = output_array[(zz * Ny + yy) * Nx + xx][1] / (Nx * Ny * Nz);
				}
			}
			MatrixOperation::matMerge(output_real, output_imag, output[zz]);
		}
	}

	else if (outputn == 1)
	{

		for (int zz = 0; zz < Nz; ++zz)
		{
			cv::Mat output_real = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);

			//std::cout << zz << std::endl;
			for (int yy = 0; yy < Ny; ++yy)
			{
				for (int xx = 0; xx < Nx; ++xx)
				{
					output_real.at<double>(yy, xx) = output_array[(zz * Ny + yy) * Nx + xx][0] / (Nx * Ny * Nz);
				}
			}

			output[zz] = output_real;
		}

	}

	fftw_free(input_array);
	fftw_free(output_array);

}