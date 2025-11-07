#include "RLD.cuh"
#include "WB_back_projector.h"

#include <omp.h>
//#include "kernel.cuh"




//const double PI = 3.141592653589793238462643383279;
//const double EPS = 1e-8;
// //选定数据计算范围
//int	Ny = 781;          
//int center_X = 1824;          
//int center_Y = 2432;           
//int Nz = 76;
//
//
//double n_m = 1;
//double lambda = 0.55;			
//double Mag = 40;				
//double pixelsize = 6.5 / Mag;	
//double z_step = 0.5;				
//double NA_obj = 0.75;			
//double k = 2 * PI / lambda;
//
////std::vector<cv::Mat> PSF(Nz), Hpsf(Nz);
//std::vector<cv::Mat> AOTF(Nz);
//std::vector<cv::Mat>Image_Ori(Nz);
//
//std::vector<cv::Mat>Object_3D_ImageEstimate(Nz);
void Matlab2C::linspace_number(double input_begin, double input_finish, int input_number, cv::Mat& Output_array)
{
	cv::Mat a1 = cv::Mat::zeros(cv::Size(1, input_number), CV_64FC1);
	double interval = (input_finish - input_begin) / (input_number - 1);
	for (int i = 0; i < input_number; i++) {
		a1.at<double>(i) = input_begin + i * interval;
	}
	Output_array = a1.clone();
}

void Matlab2C::linspace_interval(double input_begin, double input_finish, double input_interval, cv::Mat& Output_array)
{
	double number = (input_finish - input_begin) / input_interval + 1;
	cv::Mat a1 = cv::Mat::zeros(cv::Size(1, number), CV_64FC1);
	for (int i = 0; i < number; i++) {
		a1.at<double>(i) = input_begin + i * input_interval;
	}
	Output_array = a1.clone();
}

void Matlab2C::Repmat(const cv::Mat& input_mat, int repeat_num, std::vector<cv::Mat>& Output_array)
{
	cv::Mat ZerosMat;
	if (input_mat.channels() == 1) {  
		for (int i = 0; i < repeat_num; i++) {
			ZerosMat = cv::Mat::zeros(cv::Size(input_mat.rows, input_mat.cols), CV_64FC1);
			Output_array[i].push_back(ZerosMat);
		}
	}
	else if (input_mat.channels() == 2) {  
		for (int i = 0; i < repeat_num; i++) {
			ZerosMat = cv::Mat::zeros(cv::Size(input_mat.rows, input_mat.cols), CV_64FC2);
			Output_array[i].push_back(ZerosMat);
		}
	}
	for (int i = 0; i < repeat_num; i++)
	{
		Output_array[i] = input_mat;
	}
}


std::vector<cv::Mat> TransferFunction::Caculate_F(std::vector<cv::Mat>& rhoxy, std::vector<cv::Mat>& eta, double lambda,
	std::vector<cv::Mat>& rho, double rhos, double rhop)
{

	int Nz = rhoxy.size();
	int Nx = rhoxy[0].cols;
	int Ny = rhoxy[0].rows;


	std::vector<cv::Mat> m0(Nz), delta1(Nz), delta2(Nz), delta3(Nz), delta(Nz);
	for (int zz = 0; zz < Nz; zz++)
	{
		cv::divide(eta[zz], rhoxy[zz], m0[zz]);
		delta1[zz] = m0[zz].mul(eta[zz] / 2 - sqrt(1 / (lambda * lambda) - rhos * rhos));
		delta2[zz] = delta1[zz];
		delta3[zz] = m0[zz].mul(-eta[zz] / 2 - sqrt(1 / (lambda * lambda) - rhop * rhop));
		delta[zz] = cv::Mat::zeros(cv::Size(Nx, Ny), CV_64FC1);
	}

	std::vector<cv::Mat> MASK1(Nz), MASK2(Nz), MASK3(Nz);
	cv::Mat Aperture1, Aperture2, Aperture3, Aperture4, Aperture,
		m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13;


	for (int zz = 0; zz < Nz; zz++)
	{
		Aperture1 = rhoxy[zz] <= (rhop - rhos);
		Aperture1.convertTo(Aperture1, CV_64FC1, 1 / 255.0f);//double
		//MASK1[zz] = Aperture1;
		Aperture2 = rhoxy[zz] >= 0;
		Aperture2.convertTo(Aperture2, CV_64FC1, 1 / 255.0f);
		//MASK2[zz] = Aperture2;
		cv::sqrt(1 / (lambda * lambda) - (rhos - rhoxy[zz]).mul(rhos - rhoxy[zz]), m1);
		Aperture3 = (sqrt(1 / (lambda * lambda) - rhos * rhos) - m1) <= eta[zz];
		Aperture3.convertTo(Aperture3, CV_64FC1, 1 / 255.0f);
		//MASK3[zz] = Aperture3;
		cv::sqrt(1 / (lambda * lambda) - (rhos + rhoxy[zz]).mul(rhos + rhoxy[zz]), m2);
		// Value[zz] = (sqrt(1 / (lambda * lambda) - rhos * rhos) - m2);
		Aperture4 = eta[zz] <= (sqrt(1 / (lambda * lambda) - rhos * rhos) - m2);
		Aperture4.convertTo(Aperture4, CV_64FC1, 1 / 255.0f);
		//MASK4[zz] = Aperture4;
		MASK1[zz] = Aperture1.mul(Aperture2).mul(Aperture3).mul(Aperture4);
		//delta[zz] = delta1[zz].mul(MASK1[zz]);
	}

	for (int zz = 0; zz < Nz; zz++)
	{
		Aperture1 = rhoxy[zz] >= (rhop - rhos);
		Aperture1.convertTo(Aperture1, CV_64FC1, 1 / 255.0f);//double
		//MASK1[zz] = Aperture1;
		Aperture2 = rhoxy[zz] <= (rhop + rhos);
		Aperture2.convertTo(Aperture2, CV_64FC1, 1 / 255.0f);
		//MASK2[zz] = Aperture2;
		cv::sqrt(1 / (lambda * lambda) - (rhos - rhoxy[zz]).mul(rhos - rhoxy[zz]), m1);
		Aperture3 = (sqrt(1 / (lambda * lambda) - rhos * rhos) - m1) <= eta[zz];
		Aperture3.convertTo(Aperture3, CV_64FC1, 1 / 255.0f);
		//MASK3[zz] = Aperture3;
		Aperture4 = eta[zz] <= (sqrt(1 / (lambda * lambda) - rhos * rhos) - sqrt(1 / (lambda * lambda) - rhop * rhop));
		Aperture4.convertTo(Aperture4, CV_64FC1, 1 / 255.0f);
		//MASK4[zz] = Aperture4;
		MASK2[zz] = Aperture1.mul(Aperture2).mul(Aperture3).mul(Aperture4);
		//delta[zz] = delta2[zz].mul(MASK2[zz]);
	}

	for (int zz = 0; zz < Nz; zz++)
	{
		Aperture1 = rhoxy[zz] >= (rhop - rhos);
		Aperture1.convertTo(Aperture1, CV_64FC1, 1 / 255.0f);//double
		//MASK1[zz] = Aperture1;
		Aperture2 = rhoxy[zz] <= (rhop + rhos);
		Aperture2.convertTo(Aperture2, CV_64FC1, 1 / 255.0f);
		//MASK2[zz] = Aperture2;
		Aperture3 = (sqrt(1 / (lambda * lambda) - rhos * rhos) - sqrt(1 / (lambda * lambda) - rhop * rhop)) < eta[zz];
		Aperture3.convertTo(Aperture3, CV_64FC1, 1 / 255.0f);
		//MASK3[zz] = Aperture3;
		cv::sqrt(1 / (lambda * lambda) - (rhop - rhoxy[zz]).mul(rhop - rhoxy[zz]), m1);
		Aperture4 = eta[zz] <= (m1 - sqrt(1 / (lambda * lambda) - rhop * rhop));
		Aperture4.convertTo(Aperture4, CV_64FC1, 1 / 255.0f);
		//MASK4[zz] = Aperture4;
		MASK3[zz] = Aperture1.mul(Aperture2).mul(Aperture3).mul(Aperture4);
		delta[zz] = delta1[zz].mul(MASK1[zz]) + delta2[zz].mul(MASK2[zz]) + delta3[zz].mul(MASK3[zz]);
	}

	/*std::vector<cv::Mat> F(Nz), F1(Nz), mm0(Nz), mm1(Nz), mm2(Nz), mm3(Nz), mm4(Nz),
		mm5(Nz), mm6(Nz), mm7(Nz), mm8(Nz), mm9(Nz), mm10(Nz), mm11(Nz), mm12(Nz);*/

	std::vector<cv::Mat> F(Nz), F1(Nz);
	cv::Mat MASKdelta, mm0, mm1, mm2, mm3, mm4,
		mm5, mm6, mm7, mm8, mm9, mm10, mm11, mm12;

	for (int zz = 0; zz < Nz; zz++)
	{
		mm0 = rhoxy[zz].mul(rhoxy[zz]).mul(delta[zz]);
		mm1 = rho[zz].mul(rho[zz]).mul(eta[zz]);

		cv::divide(mm0, mm1, mm2);
		mm3 = cv::abs(mm2); // m3 p1 abs

		mm4 = 1 / (lambda * lambda) - rho[zz].mul(rho[zz]) / 4;
		cv::divide(rho[zz].mul(rho[zz]).mul(delta[zz]).mul(delta[zz]), (eta[zz].mul(eta[zz])), mm5);
		cv::sqrt(mm4 - mm5, mm6); // m6  p2 sqrt
	cv:divide(eta[zz].mul(eta[zz]), 2 * rho[zz], mm7);
		mm8 = mm4 - mm7;  // m8 p3
		cv::sqrt(mm4, mm9);
		cv::divide(rho[zz].mul(delta[zz]), eta[zz].mul(mm9), mm10);
		mm11 = cv::abs(mm10);
		MatrixOperation::Mat_acos(mm11, mm12);
		F1[zz].push_back(mm3.mul(mm6) + mm8.mul(mm12));
		MASKdelta = ~(delta[zz] == 0);
		MASKdelta.convertTo(MASKdelta, CV_64FC1, 1 / 255.0f);//double
		F[zz] = F1[zz].mul(MASKdelta);

		//mm0[zz] = rhoxy[zz].mul(rhoxy[zz]).mul(delta[zz]);
		//mm1[zz] = rho[zz].mul(rho[zz]).mul(eta[zz]);

		//cv::divide(mm0[zz], mm1[zz], mm2[zz]);
		//mm3[zz] = cv::abs(mm2[zz]); // m3 p1 abs

		//mm4[zz] = 1 / (lambda * lambda) - rho[zz].mul(rho[zz]) / 4;
		//cv::divide(rho[zz].mul(rho[zz]).mul(delta[zz]).mul(delta[zz]), (eta[zz].mul(eta[zz])), mm5[zz]);
		//cv::sqrt(mm4[zz] - mm5[zz], mm6[zz]); // m6  p2 sqrt
		//cv:divide(eta[zz].mul(eta[zz]), 2 * rho[zz], mm7[zz]);
		//mm8[zz] = mm4[zz] - mm7[zz];  // m8 p3
		//cv::sqrt(mm4[zz], mm9[zz]);
		//cv::divide(rho[zz].mul(delta[zz]), eta[zz].mul(mm9[zz]), mm10[zz]);
		//mm11[zz] = cv::abs(mm10[zz]);
		//MatrixOperation::Mat_acos(mm11[zz], mm12[zz]);
		//F1[zz] = mm3[zz].mul(mm6[zz]) + mm8[zz].mul(mm12[zz]);
		//MASKdelta = ~(delta[zz] == 0);
		//MASKdelta.convertTo(MASKdelta, CV_64FC1, 1 / 255.0f);//double
		//F[zz] = F1[zz].mul(MASKdelta);

	}
	return F;

}



void MatrixOperation::meshgrid(cv::Mat input_x, cv::Mat input_y, cv::Mat& X, cv::Mat& Y)  
{
	cv::repeat(input_x.t(), input_x.rows, 1, X); 
	cv::repeat(input_y, 1, input_y.rows, Y);
}

void MatrixOperation::exp_ix(const cv::Mat& input_x, cv::Mat& Output_real, cv::Mat& Output_imag) { 
	cv::Mat m1;
	m1 = cv::Mat::ones(cv::Size(input_x.rows, input_x.cols), CV_64FC1);
	cv::polarToCart(m1, input_x, Output_real, Output_imag);//Output_real->cos(x);  Output_imag->sin(x);
}

void MatrixOperation::circshift(cv::Mat& input_array, const cv::Point& move_num)
{
	cv::Size sz = input_array.size();

	
	assert(sz.height > 0 && sz.width > 0); 

	if ((sz.height == 1 && sz.width == 1) || (move_num.x == 0 && move_num.y == 0))
		return;


	int x = move_num.x;
	int y = move_num.y;
	if (x > 0) x = x % sz.width;
	if (y > 0) y = y % sz.height;
	if (x < 0) x = x % sz.width + sz.width;
	if (y < 0) y = y % sz.height + sz.height;


	std::vector<cv::Mat> planes;
	split(input_array, planes);

	for (size_t i = 0; i < planes.size(); i++)
	{
		cv::Mat tmp0, tmp1, tmp2, tmp3;

		if (y != 0) {
			cv::Mat q0(planes[i], cv::Rect(0, 0, sz.width, sz.height - y)); //Rect(INT, INT, INT, INT)
			cv::Mat q1(planes[i], cv::Rect(0, sz.height - y, sz.width, y));
			q0.copyTo(tmp0);
			q1.copyTo(tmp1);
			tmp0.copyTo(planes[i](cv::Rect(0, y, sz.width, sz.height - y)));
			tmp1.copyTo(planes[i](cv::Rect(0, 0, sz.width, y)));
		}


		if (x != 0) {
			cv::Mat q2(planes[i], cv::Rect(0, 0, sz.width - x, sz.height));
			cv::Mat q3(planes[i], cv::Rect(sz.width - x, 0, x, sz.height));
			q2.copyTo(tmp2);
			q3.copyTo(tmp3);
			tmp2.copyTo(planes[i](cv::Rect(x, 0, sz.width - x, sz.height)));
			tmp3.copyTo(planes[i](cv::Rect(0, 0, x, sz.height)));
		}

	}
	merge(planes, input_array);
}

void MatrixOperation::matSplit(const cv::Mat& input_array, cv::Mat& Output_real, cv::Mat& Output_img) 
{
	cv::Mat planes[2];
	cv::split(input_array, planes);
	(planes[0]).copyTo(Output_real);//
	(planes[1]).copyTo(Output_img);//

}

void MatrixOperation::matMerge(const cv::Mat& input_real, const cv::Mat& input_imag, cv::Mat& Output_merge)
{
	cv::Mat planes[2];
	input_real.copyTo(planes[0]);
	input_imag.copyTo(planes[1]);
	cv::merge(planes, 2, Output_merge);
}


void MatrixOperation::cMatMultiply(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array) 
{

	int row1 = input_array1.rows;//
	int col1 = input_array1.cols;

	int row2 = input_array2.rows;//
	int col2 = input_array2.cols;

	int row_d = Output_array.rows;//
	int col_d = Output_array.cols;


	assert(col1 == col2 && row1 == row2); //assert 

	
	cv::Mat src1_real, src1_imag, src2_real, src2_imag;
	cv::Mat channel[2] = { cv::Mat::zeros(cv::Size(col_d,row_d), CV_64FC1) ,cv::Mat::zeros(cv::Size(col_d, row_d), CV_64FC1) };

	matSplit(input_array1, src1_real, src1_imag);
	matSplit(input_array2, src2_real, src2_imag);


	channel[0] = src1_real.mul(src2_real) - src1_imag.mul(src2_imag);
	channel[1] = src1_imag.mul(src2_real) + src1_real.mul(src2_imag);


	merge(channel, 2, Output_array);
}

void MatrixOperation::Mat_Mul_cMat(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array)

{
	cv::Mat input_array2_real, input_array2_imag;
	matSplit(input_array2, input_array2_real, input_array2_imag);
	cv::Mat Output_array_real, Output_array_imag;

	Output_array_real = input_array1.mul(input_array2_real);
	Output_array_imag = input_array1.mul(input_array2_imag);
	matMerge(Output_array_real, Output_array_imag, Output_array);
}


void MatrixOperation::cMatDivide(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array)
{
	
	int row1 = input_array1.rows;
	int col1 = input_array1.cols;

	int row2 = input_array2.rows;
	int col2 = input_array2.cols;

	int row_d = Output_array.rows;
	int col_d = Output_array.cols;


	assert(col1 == col2 && row1 == row2); 


	cv::Mat src1_real, src1_imag, src2_real, src2_imag;
	cv::Mat channel[2] = { cv::Mat::zeros(cv::Size(col_d,row_d), CV_64FC1) ,cv::Mat::zeros(cv::Size(col_d, row_d), CV_64FC1) };

	matSplit(input_array1, src1_real, src1_imag);//a+bi
	matSplit(input_array2, src2_real, src2_imag);//c+di

	// (a+bi)./(c+di)=(ac+bd/(c.^2+d.^2))+(bc-ad/(c.^2+d.^2))i = ac+bd/(c.^2+d.^2) + bc-ad/(c.^2+d.^2)i
	channel[0] = (src1_real.mul(src2_real) + src1_imag.mul(src2_imag)) / (src2_real.mul(src2_real) + src2_imag.mul(src2_imag));
	channel[1] = (src1_imag.mul(src2_real) - src1_real.mul(src2_imag)) / (src2_real.mul(src2_real) + src2_imag.mul(src2_imag));

	
	merge(channel, 2, Output_array);
}

void MatrixOperation::cMat_Div_Mat(const cv::Mat& input_array1, const cv::Mat& input_array2, cv::Mat& Output_array)

{
	cv::Mat input_array2_real, input_array2_imag;
	matSplit(input_array2, input_array2_real, input_array2_imag);
	cv::Mat Output_array_real, Output_array_imag;

	Output_array_real = input_array1 / input_array2_real;
	Output_array_imag = input_array1 / input_array2_imag;
	matMerge(Output_array_real, Output_array_imag, Output_array);
}

void MatrixOperation::cMat_Div_Mat3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output)
{
	int Nz = input1.size();
	int Nx = input1[0].rows;
	int Ny = input1[0].cols;
	int channel = input1[0].channels();

	if (channel == 1)
	{
		cv::Mat input2_real, input2_imag, output_real, output_imag;

		for (int zz = 0; zz < Nz; zz++) {

			matSplit(input2[zz], input2_real, input2_imag);
			output_real = input1[zz] / (input2_real);
			output_imag = input1[zz] / (input2_imag);
			matMerge(output_real, output_imag, output[zz]);

		}
	}
	else if (channel == 2)
	{
		cv::Mat input1_real, input1_imag, input2_real, input2_imag, output_real, output_imag;

		for (int zz = 0; zz < Nz; zz++) {

			matSplit(input1[zz], input1_real, input1_imag);
			matSplit(input2[zz], input2_real, input2_imag);
			//(a+bi)./(c+di)=(ac+bd/(c.^2+d.^2))+(bc-ad/(c.^2+d.^2))i = (ac+bd) / (c.^2+d.^2) + (bc-ad) / (c.^2+d.^2)i
			output_real = (input1_real.mul(input2_real) + input1_imag.mul(input2_imag)) / (input2_real.mul(input2_real) + input2_imag.mul(input2_imag));
			output_imag = (input1_imag.mul(input2_real) - input1_real.mul(input2_imag)) / (input2_real.mul(input2_real) + input2_imag.mul(input2_imag));

			matMerge(output_real, output_imag, output[zz]);

		}

	}


}



void MatrixOperation::cMat_Div3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output)
{
	int Nz = input1.size();



	for (int zz = 0; zz < Nz; zz++) {

		output[zz] = input1[zz] / input2[zz];

	}


}


void MatrixOperation::cMatAdd(cv::Mat& input_array1, cv::Mat& input_array2, cv::Mat& Output_array)//
{

	int row1 = input_array1.rows;
	int col1 = input_array1.cols;

	int row2 = input_array2.rows;
	int col2 = input_array2.cols;

	int row_d = Output_array.rows;
	int col_d = Output_array.cols;

	
	assert(col1 == col2 && row1 == row2); //assert 


	cv::Mat src1_real, src1_imag, src2_real, src2_imag;
	cv::Mat channel[2] = { cv::Mat::zeros(cv::Size(col_d,row_d), CV_64FC1) ,cv::Mat::zeros(cv::Size(col_d, row_d), CV_64FC1) };

	matSplit(input_array1, src1_real, src1_imag);
	matSplit(input_array2, src2_real, src2_imag);

	// (a+bi)+(c+di)=（a+c)+(b+d)i = a+c + ：b+d
	channel[0] = src1_real + src2_real;
	channel[1] = src1_imag + src2_imag;

	
	merge(channel, 2, Output_array);
}


void MatrixOperation::Mat_acos(cv::Mat& input_array, cv::Mat& Output_array)
{

	int w = input_array.cols;
	int h = input_array.rows;
	Output_array = cv::Mat::zeros(cv::Size(w, h), CV_64FC1);
	for (int row = 0; row < h; row++)
	{
		for (int col = 0; col < w; col++)
		{
			double pv = input_array.at<double>(row, col); // 

			Output_array.at<double>(row, col) = acos(pv); // 

		}
	}



}


void MatrixOperation::mapminmax(const cv::Mat& input_array, double max_Value, double min_Value, cv::Mat& Output_array)
{
	double input_array_min = 0, input_array_max = 0;
	cv::minMaxLoc(input_array, &input_array_min, &input_array_max, 0, 0);
	cv::Mat y;
	y = (max_Value - min_Value) * (input_array - input_array_min) / (input_array_max - input_array_min) + min_Value;
	y.copyTo(Output_array);
}

void MatrixOperation::Nonnegative(cv::Mat& input_array)
{
	int w = input_array.cols;
	int h = input_array.rows;

	for (int row = 0; row < h; row++)
	{
		for (int col = 0; col < w; col++)
		{
			double pv = input_array.at<double>(row, col); // 
			if (pv < 0) {
				input_array.at<double>(row, col) = 0; // 
			}

		}
	}



}


void MatrixOperation::Mat_Mul3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output)
{
	int Nz = input1.size();
	for (int zz = 0; zz < Nz; zz++) {

		output[zz] = input1[zz].mul(input2[zz]);
	}
}

void MatrixOperation::Mat_Mul_cMat3(const std::vector<cv::Mat>& input1, const std::vector<cv::Mat>& input2, std::vector<cv::Mat>& output)
{
	int Nz = input1.size();
	int Nx = input1[0].rows;
	int Ny = input1[0].cols;
	int channel = input1[0].channels();

	if (channel == 1)
	{
		cv::Mat input2_real, input2_imag, output_real, output_imag;

		for (int zz = 0; zz < Nz; zz++) {

			matSplit(input2[zz], input2_real, input2_imag);
			output_real = input1[zz].mul(input2_real);
			output_imag = input1[zz].mul(input2_imag);
			matMerge(output_real, output_imag, output[zz]);

		}
	}
	else if (channel == 2)
	{
		cv::Mat input1_real, input1_imag, input2_real, input2_imag, output_real, output_imag;

		for (int zz = 0; zz < Nz; zz++) {

			matSplit(input1[zz], input1_real, input1_imag);
			matSplit(input2[zz], input2_real, input2_imag);
			output_real = input1_real.mul(input2_real) - input1_imag.mul(input2_imag);
			output_imag = input1_imag.mul(input2_real) + input1_real.mul(input2_imag);

			matMerge(output_real, output_imag, output[zz]);

		}
	}
}

__global__ void kernel_copy(cufftComplex* data1, cufftComplex* data2, int data_len)//real_mul data2 output
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {

		data2[idx].x = data1[idx].x;
		data2[idx].y = data1[idx].y;
	}

}
__global__ void kernel_copy1(cufftComplex* data1, float* data2, int data_len)//real_mul data2 output
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {

		data2[idx] = data1[idx].x;
		//data2[idx].y = data1[idx].y;
	}

}
__global__ void kernel_cMat_Div3(cufftComplex* data1, cufftComplex* data2, int data_len)//real div
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {
		//data1[idx].x	input1.real	// data1[idx].y	input1.imag
		//data2[idx].x	input2.real	//data2[idx].y	input2.imag
		data2[idx].x = data1[idx].x / data2[idx].x;

	}

}
__global__ void kernel_Mat_Mul3(cufftComplex* data1, cufftComplex* data2, int data_len)//real_mul data2 output
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {

		data2[idx].x = data1[idx].x * data2[idx].x;

	}

}
__global__ void kernel_Mat_Mul_cMat3(cufftComplex* data1, cufftComplex* data2, int data_len)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {
		//data1[idx].x	input1.real	// data1[idx].y	input1.imag
		//data2[idx].x	input2.real	//data2[idx].y	input2.imag
		data2[idx].x = data1[idx].x * data2[idx].x - data1[idx].y * data2[idx].y;
		data2[idx].y = data1[idx].x * data2[idx].y + data1[idx].y * data2[idx].x;
	}

}
__global__ void kernel_Nonnegative(cufftComplex* data, int data_len)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {
		if (data[idx].x < 0)
		{
			data[idx].x = (float)0;
		}
		if (data[idx].y < 0)
		{
			data[idx].y = (float)0;
		}

	}

}
__global__ void kernel_normalizing(cufftComplex* data, int data_len)//
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {
		data[idx].x /= (float)data_len;
		data[idx].y /= (float)data_len;
	}

}
__global__ void kernel_normalizing_imag0(cufftComplex* data, int data_len)//
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {
		data[idx].x /= (float)data_len;
		data[idx].y = (float)0;
	}

}
__global__ void kernel_imag0(cufftComplex* data, int data_len)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < data_len) {

		data[idx].y = (float)0;
	}

}
__global__ void kernel_fftshift3(cufftComplex* data, cufftComplex* out, int X, int Y, int Z)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int Size = X * Y * Z;
	if (idx < Size)
	{
		int z = idx % Z;
		int xy = idx / Z;
		int y = xy % Y;
		int x = xy / Y;

		
		int shiftX = (x + X / 2) % X;
		int shiftY = (y + Y / 2) % Y;
		int shiftZ = (z + Z / 2) % Z;
		
		int shiftIdx = ((shiftX * Y + shiftY) * Z + shiftZ);
		out[shiftIdx] = data[idx];

	}
}
__global__ void kernel_ifftshift3(cufftComplex* data, cufftComplex* out, int X, int Y, int Z)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int Size = X * Y * Z;
	if (idx < Size)
	{
		int z = idx % Z;
		int xy = idx / Z;
		int y = xy % Y;
		int x = xy / Y;

		
		int shiftX = (x + (X + 1) / 2) % X;
		int shiftY = (y + (Y + 1) / 2) % Y;
		int shiftZ = (z + (Z + 1) / 2) % Z;
		
		int shiftIdx = ((shiftX * Y + shiftY) * Z + shiftZ);
		out[shiftIdx] = data[idx];

	}
}
void normalizeImages(std::vector<cv::Mat>& images) {
	for (auto& image : images) {
		cv::normalize(image, image, 0, 1, cv::NORM_MINMAX);
	}
}
void globalNormalize(std::vector<cv::Mat>& images) {
	
	double globalMin = DBL_MAX;
	double globalMax = DBL_MIN;

	for (const auto& image : images) {
		double minVal, maxVal;
		cv::minMaxLoc(image, &minVal, &maxVal);
		globalMin = std::min(globalMin, minVal);
		globalMax = std::max(globalMax, maxVal);
	}

	
	for (auto& image : images) {
		image.convertTo(image, CV_64F); 
		image = (image - globalMin) / (globalMax - globalMin);
	}
}




void parallelDataTransfer(const std::vector<cv::Mat>& source,
	cufftComplex* dest,
	int Nx, int Ny, int Nz)
{
#pragma omp parallel for collapse(2)
	for (int k = 0; k < Nz; k++) {
		for (int i = 0; i < Nx; i++) {
			cv::Mat& non_const_mat = const_cast<cv::Mat&>(source[k]);
			double* row_ptr = non_const_mat.ptr<double>(i);
			for (int j = 0; j < Ny; j++) {
				size_t index = i * Ny * Nz + j * Nz + k;
				dest[index].x = row_ptr[j];
				dest[index].y = 0;
			}
		}
	}
}


void copyFromGPU(cufftComplex* d_Object_3D,
	std::vector<cv::Mat>& Object_3D_ImageEstimate,
	int Nx, int Ny, int Nz)
{
	
	Object_3D_ImageEstimate.resize(Nz);

#pragma omp parallel for
	for (int k = 0; k < Nz; k++) {
		cv::Mat mat(Nx, Ny, CV_64FC1);

		
		std::vector<double*> row_ptrs(Nx);
		for (int i = 0; i < Nx; i++) {
			row_ptrs[i] = mat.ptr<double>(i);
		}

		
#pragma omp simd collapse(2)
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				size_t index = i * Ny * Nz + j * Nz + k;
				row_ptrs[i][j] = d_Object_3D[index].x;
			}
		}

	
		Object_3D_ImageEstimate[k] = mat;
	}
}




int NativeRLDeconvolution(const std::vector<cv::Mat>&Image_Ori,const std::vector<cv::Mat>&AOTF,  std::vector<cv::Mat>& Object_3D_ImageEstimate , int NI)//RLD
{
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	

	cudaEventRecord(start);

	int Nx = cv::Mat(Image_Ori[0]).rows;
	int Ny = cv::Mat(Image_Ori[0]).cols;
	int Nz = Image_Ori.size();

	cufftHandle plan;
	cufftComplex* idata, * odata;
	cufftComplex* idata1, * odata1;
	cufftComplex* ImageEstimate, * d_Blurred_3D;
	//float* dst;
	cudaMalloc(&ImageEstimate, sizeof(cufftComplex) * Nx * Ny * Nz);//
	cudaMalloc(&d_Blurred_3D, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&idata, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&odata, sizeof(cufftComplex) * Nx * Ny * Nz);
	//cudaMalloc(&dst, sizeof(float)* Nx* Ny* Nz);

	cudaMalloc(&idata1, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&odata1, sizeof(cufftComplex) * Nx * Ny * Nz);

	cufftComplex* d_AOTF = (cufftComplex*)malloc(sizeof(cufftComplex) * Nx * Ny * Nz);//AOTF
	cufftComplex* d_Object_3D = (cufftComplex*)malloc(sizeof(cufftComplex) * Nx * Ny * Nz);//Object_3D
	int n[3] = { Nx, Ny, Nz };
	cufftPlanMany(&plan, 3, n,
		NULL, 1, Nx * Ny * Nz,
		NULL, 1, Nx * Ny * Nz,
		CUFFT_C2C, 1);

	for (int k = 0; k < Nz; k++) {//AOTF
		cv::Mat mat = AOTF[k];
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				cv::Vec2f vec = mat.at<cv::Vec2f>(i, j);
				d_AOTF[i * Ny * Nz + j * Nz + k].x = AOTF[k].at<double>(i, j);
				d_AOTF[i * Ny * Nz + j * Nz + k].y = 0;
				/*	cufftComplex vec = mat.at<cufftComplex>(i, j);
					d_AOTF[i * Ny * Nz + j * Nz + k] = vec;*/
			}
		}
	}
	for (int k = 0; k < Nz; k++) {//Object_3D
		cv::Mat mat = Image_Ori[k];
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				cv::Vec2f vec = mat.at<cv::Vec2f>(i, j);
				d_Object_3D[i * Ny * Nz + j * Nz + k].x = Image_Ori[k].at<double>(i, j);//Z_array		Object_3D
				d_Object_3D[i * Ny * Nz + j * Nz + k].y = 0;
				/*	cufftComplex vec = mat.at<cufftComplex>(i, j);
					d_AOTF[i * Ny * Nz + j * Nz + k] = vec;*/
			}
		}
	}


	/////////////
	//cudaMalloc((void**)&idata, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMemcpy(idata, d_AOTF, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyHostToDevice);//idata (AOTF) host device
	cudaMemcpy(idata1, d_Object_3D, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyHostToDevice);//idata1 (d_Object_3D)

	//block	grid
	int block_size = 256;
	int totalSize = Nx * Ny * Nz;
	int grid_size = (totalSize + block_size - 1) / block_size;

	kernel_ifftshift3 <<<grid_size, block_size >>> (idata, odata, Nx, Ny, Nz);
	cufftExecC2C(plan, odata, odata, CUFFT_INVERSE);//CUFFT_FORWARD
	kernel_normalizing_imag0 << <grid_size, block_size >> > (odata, totalSize); //PSF			
	kernel_fftshift3 << <grid_size, block_size >> > (odata, idata, Nx, Ny, Nz);//PSFa	

	kernel_Nonnegative << <grid_size, block_size >> > (idata, totalSize); //	

	kernel_ifftshift3 << <grid_size, block_size >> > (idata, odata, Nx, Ny, Nz);//			
	cufftExecC2C(plan, odata, odata, CUFFT_FORWARD);//Hpsf	

	//kernel_normalizing_imag0 << <grid_size, block_size >> > (odata, totalSize);
	kernel_fftshift3 << <grid_size, block_size >> > (odata, idata, Nx, Ny, Nz);	//				
	kernel_imag0 << <grid_size, block_size >> > (idata, totalSize);//idata	Hpsf			
	//////idata保存		Hpsf	1


	ImageEstimate = idata1;

	kernel_copy << <grid_size, block_size >> > (idata1, d_Blurred_3D, totalSize);
	//cudaMemcpy(d_AOTF, idata1, sizeof(cufftComplex)* Nx* Ny* Nz, cudaMemcpyDeviceToHost);
	/*int NI = 200;*/
	//迭代
	for (int niter = 0; niter < NI; niter++)
	{
		std::cout << "niter = " << niter << std::endl;
		cufftExecC2C(plan, ImageEstimate, odata1, CUFFT_FORWARD);//HI	odata1									

		kernel_fftshift3 << <grid_size, block_size >> > (odata1, odata, Nx, Ny, Nz);//HI	odata			

		kernel_Mat_Mul_cMat3 << <grid_size, block_size >> > (idata, odata, totalSize);//HConv	odata				

		kernel_ifftshift3 << <grid_size, block_size >> > (odata, odata1, Nx, Ny, Nz);//HConv	odata1			

		cufftExecC2C(plan, odata1, odata1, CUFFT_INVERSE);		//Conv	odata1									

		kernel_normalizing_imag0 << <grid_size, block_size >> > (odata1, totalSize);//Conv	odata1		

		kernel_cMat_Div3 << <grid_size, block_size >> > (d_Blurred_3D, odata1, totalSize);//DV	odata1	 right		

		cufftExecC2C(plan, odata1, odata1, CUFFT_FORWARD);//HDV		odata1										

		kernel_fftshift3 << <grid_size, block_size >> > (odata1, odata, Nx, Ny, Nz);//HDV odata				
		//
		kernel_Mat_Mul_cMat3 << <grid_size, block_size >> > (idata, odata, totalSize);//HDV_Conv	odata		

		kernel_ifftshift3 << <grid_size, block_size >> > (odata, odata1, Nx, Ny, Nz);//	HDV_Conv	odata1		

		cufftExecC2C(plan, odata1, odata1, CUFFT_INVERSE);//

		kernel_normalizing_imag0 << <grid_size, block_size >> > (odata1, totalSize);//DV_Conv	odata1	1		

		kernel_Mat_Mul3 << <grid_size, block_size >> > (odata1, ImageEstimate, totalSize);
		cudaDeviceSynchronize();
	}
	//kernel_copy1 << <grid_size, block_size >> > (ImageEstimate, dst, totalSize);

	//
	cudaMemcpy(d_Object_3D, ImageEstimate, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	cufftDestroy(plan);
	cudaFree(idata);
	cudaFree(odata);
	cudaFree(idata1);
	cudaFree(odata1);

	for (int k = 0; k < Nz; k++) {
		cv::Mat mat(Nx, Ny, CV_64FC1);
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {

				mat.at< double>(i, j) = d_Object_3D[i * Ny * Nz + j * Nz + k].x;

			}
		}
		Object_3D_ImageEstimate[k] = mat;//Object_3D
	}
	
	cudaEventRecord(stop);
	cudaEventSynchronize(stop); //
	
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	
	std::cout << "迭代: " << NI << " RLD_cuda总计算时间: " << milliseconds << " ms" << std::endl;
	std::cout << "平均每次迭代时间，含数据传输: " << milliseconds / NI << " ms" << std::endl;
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	normalizeImages(Object_3D_ImageEstimate);

	cv::Mat ImageEstimate1 = Object_3D_ImageEstimate[0];//
	cv::Mat ImageEstimate15 = Object_3D_ImageEstimate[14];//
	cv::Mat ImageEstimate26 = Object_3D_ImageEstimate[25];//
	cv::Mat ImageEstimate27 = Object_3D_ImageEstimate[26];//
	cv::Mat ImageEstimate41 = Object_3D_ImageEstimate[40];//
	cv::Mat ImageEstimate45 = Object_3D_ImageEstimate[44];//
	cv::Mat ImageEstimate50 = Object_3D_ImageEstimate[49];//
	cv::Mat ImageEstimate60 = Object_3D_ImageEstimate[59];//



	for (int i = 0; i < Nz; i++)
	{
		cv::imshow("Object_3D", Image_Ori[i]);
		cv::imshow("ImageEstimate", Object_3D_ImageEstimate[i]);
		cv::waitKey(100);
	}

	//std::string path = "D://openCV//opencv_file//3DRLD_TEST//RLD//3D RLD//result//";
	//for (int i = 0; i < Nz; i++) {
	//	std::string filename;

	//	cv::Mat img = Object_3D_ImageEstimate[i];
	//	cv::Mat img_save;
	//	img.convertTo(img_save, CV_32FC1);
	//	filename = path + std::to_string(i + 1) + ".tif";
	//	cv::imwrite(filename, img_save);
	//}
	return 0;
}


int WienerButterworth_ARLDeconvolution( std::vector<cv::Mat>& Image_Ori,  std::vector<cv::Mat>& PSF_fp,  std::vector<cv::Mat>& PSF_bp, std::vector<cv::Mat>& Object_3D_ImageEstimate , int NI)//WBRLD
{
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);


	cudaEventRecord(start);
	int Nx = cv::Mat(Image_Ori[0]).rows;
	int Ny = cv::Mat(Image_Ori[0]).cols;
	int Nz = Image_Ori.size();

	cufftHandle plan;
	cufftComplex* idata, * odata;
	cufftComplex* idata_PSF_bp;
	cufftComplex* idata1, * odata1;
	cufftComplex* ImageEstimate, * d_Blurred_3D;
	
	//float* dst;
	cudaMalloc(&ImageEstimate, sizeof(cufftComplex) * Nx * Ny * Nz);//
	cudaMalloc(&d_Blurred_3D, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&idata, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&odata, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&idata_PSF_bp, sizeof(cufftComplex) * Nx * Ny * Nz);
	//cudaMalloc(&dst, sizeof(float)* Nx* Ny* Nz);

	cudaMalloc(&idata1, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMalloc(&odata1, sizeof(cufftComplex) * Nx * Ny * Nz);
	
	cufftComplex* d_PSF_fp = (cufftComplex*)malloc(sizeof(cufftComplex) * Nx * Ny * Nz);//PSF_fp
	cufftComplex* d_PSF_bp = (cufftComplex*)malloc(sizeof(cufftComplex) * Nx * Ny * Nz);//PSF_bp
	cufftComplex* d_Object_3D = (cufftComplex*)malloc(sizeof(cufftComplex) * Nx * Ny * Nz);//Object_3D
	int n[3] = { Nx, Ny, Nz };
	cufftPlanMany(&plan, 3, n,
		NULL, 1, Nx * Ny * Nz,
		NULL, 1, Nx * Ny * Nz,
		CUFFT_C2C, 1);


	
	
	omp_set_num_threads(omp_get_max_threads());

	//omp_set_dynamic(0);                      
	//omp_set_num_threads(24);               


	// PSF_fp
	parallelDataTransfer(PSF_fp, d_PSF_fp, Nx, Ny, Nz);

	// PSF_bp
	parallelDataTransfer(PSF_bp, d_PSF_bp, Nx, Ny, Nz);

	// Image_Ori
	parallelDataTransfer(Image_Ori, d_Object_3D, Nx, Ny, Nz);



	//cudaMalloc((void**)&idata, sizeof(cufftComplex) * Nx * Ny * Nz);
	cudaMemcpy(idata, d_PSF_fp, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyHostToDevice);//idata (PSF_fp) host device
	cudaMemcpy(idata1, d_Object_3D, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyHostToDevice);//idata1 (d_Object_3D)
	cudaMemcpy(idata_PSF_bp, d_PSF_bp, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyHostToDevice);//idata_PSF_bp (d_PSF_bp)
	//block	grid
	int block_size = 256;
	int totalSize = Nx * Ny * Nz;
	int grid_size = (totalSize + block_size - 1) / block_size;

	//kernel_ifftshift3 << <grid_size, block_size >> > (idata, odata, Nx, Ny, Nz);
	//cufftExecC2C(plan, odata, odata, CUFFT_INVERSE);//CUFFT_FORWARD
	//kernel_normalizing_imag0 << <grid_size, block_size >> > (odata, totalSize); //PSF			
	//kernel_fftshift3 << <grid_size, block_size >> > (odata, idata, Nx, Ny, Nz);//PSFa	

	//kernel_Nonnegative << <grid_size, block_size >> > (idata, totalSize); //PSFa 

	kernel_ifftshift3 << <grid_size, block_size >> > (idata, odata, Nx, Ny, Nz);//			
	cufftExecC2C(plan, odata, odata, CUFFT_FORWARD);//Hpsf	
	//kernel_normalizing_imag0 << <grid_size, block_size >> > (odata, totalSize);
	kernel_fftshift3 << <grid_size, block_size >> > (odata, idata, Nx, Ny, Nz);	//	
	kernel_imag0 << <grid_size, block_size >> > (idata, totalSize);//idata	Hpsf			
	//////idata	Hpsf	1

	/////////////////////////////     odata1
	kernel_ifftshift3 << <grid_size, block_size >> > (idata_PSF_bp, odata1, Nx, Ny, Nz);//idata_PSF_bp	odata1			ifftshift	
	cufftExecC2C(plan, odata1, odata1, CUFFT_FORWARD);//HPSFpFlip	 fftn
	//kernel_normalizing_imag0 << <grid_size, block_size >> > (odata, totalSize);
	kernel_fftshift3 << <grid_size, block_size >> > (odata1, idata_PSF_bp, Nx, Ny, Nz);	//	fftshift
	kernel_imag0 << <grid_size, block_size >> > (idata_PSF_bp, totalSize);//HPSFpFlip		idata_PSF_bp				
	//////idata保存		Hpsf	1


	ImageEstimate = idata1;

	kernel_copy << <grid_size, block_size >> > (idata1, d_Blurred_3D, totalSize);
	//cudaMemcpy(d_AOTF, idata1, sizeof(cufftComplex)* Nx* Ny* Nz, cudaMemcpyDeviceToHost);



	/*int NI = 2;*/
	//
	for (int niter = 0; niter < NI; niter++)
	{
		std::cout << "niter = " << niter << std::endl;
		cufftExecC2C(plan, ImageEstimate, odata1, CUFFT_FORWARD);//HI	odata1									

		kernel_fftshift3 << <grid_size, block_size >> > (odata1, odata, Nx, Ny, Nz);//HI	odata


		kernel_Mat_Mul_cMat3 << <grid_size, block_size >> > (idata, odata, totalSize);//HConv	odata			idata是Hpsf 	 odata是HI			odata是HConv

		kernel_ifftshift3 << <grid_size, block_size >> > (odata, odata1, Nx, Ny, Nz);//HConv	odata1			

		cufftExecC2C(plan, odata1, odata1, CUFFT_INVERSE);		//Conv	odata1									

		kernel_normalizing_imag0 << <grid_size, block_size >> > (odata1, totalSize);//Conv	odata1


		kernel_cMat_Div3 << <grid_size, block_size >> > (d_Blurred_3D, odata1, totalSize);//DV	odata1	 right	


		cufftExecC2C(plan, odata1, odata1, CUFFT_FORWARD);//HDV		odata1										
		kernel_fftshift3 << <grid_size, block_size >> > (odata1, odata, Nx, Ny, Nz);//HDV odata	

		//data HPSFpFlip
		kernel_Mat_Mul_cMat3 << <grid_size, block_size >> > (idata_PSF_bp, odata, totalSize);//HDV_Conv	odata		

		kernel_ifftshift3 << <grid_size, block_size >> > (odata, odata1, Nx, Ny, Nz);//	HDV_Conv	odata1		

		cufftExecC2C(plan, odata1, odata1, CUFFT_INVERSE);//

		kernel_normalizing_imag0 << <grid_size, block_size >> > (odata1, totalSize);//DV_Conv	odata1	1		

		kernel_Mat_Mul3 << <grid_size, block_size >> > (odata1, ImageEstimate, totalSize);
		kernel_Nonnegative << <grid_size, block_size >> > (ImageEstimate, totalSize);
		cudaDeviceSynchronize();
	}
	//kernel_copy1 << <grid_size, block_size >> > (ImageEstimate, dst, totalSize);



	//
	cudaMemcpy(d_Object_3D, ImageEstimate, sizeof(cufftComplex) * Nx * Ny * Nz, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	cufftDestroy(plan);
	cudaFree(idata);
	cudaFree(odata);
	cudaFree(idata1);
	cudaFree(odata1);



	//for (int k = 0; k < Nz; k++) {  ///////  
	//	cv::Mat mat(Nx, Ny, CV_64FC1);
	//	for (int i = 0; i < Nx; i++) {
	//		for (int j = 0; j < Ny; j++) {

	//			mat.at< double>(i, j) = d_Object_3D[i * Ny * Nz + j * Nz + k].x;

	//		}
	//	}
	//	Object_3D_ImageEstimate[k] = mat;//Object_3D
	//}
	
	copyFromGPU(d_Object_3D, Object_3D_ImageEstimate, Nx, Ny, Nz);

	// 
	cudaEventRecord(stop);
	cudaEventSynchronize(stop); //
	// 
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	// 
	std::cout << "迭代: " << NI << " WB-ARL_cuda总计算时间: " << milliseconds << " ms" << std::endl;
	std::cout << "平均每次迭代时间，含数据传输: " << milliseconds / NI << " ms" << std::endl;
	// 
	cudaEventDestroy(start);
	cudaEventDestroy(stop);


	//normalizeImages(Object_3D_ImageEstimate);
	//globalNormalize(Object_3D_ImageEstimate);
	cv::Mat ImageEstimate1 = Object_3D_ImageEstimate[0];//
	cv::Mat ImageEstimate15 = Object_3D_ImageEstimate[14];//
	cv::Mat ImageEstimate26 = Object_3D_ImageEstimate[25];//
	cv::Mat ImageEstimate27 = Object_3D_ImageEstimate[26];//
	cv::Mat ImageEstimate41 = Object_3D_ImageEstimate[40];//
	cv::Mat ImageEstimate45 = Object_3D_ImageEstimate[44];//
	cv::Mat ImageEstimate50 = Object_3D_ImageEstimate[49];//
	cv::Mat ImageEstimate60 = Object_3D_ImageEstimate[59];//



	for (int i = 0; i < Nz; i++)
	{
		cv::imshow("Object_3D", Image_Ori[i]);
		cv::imshow("WB-ARL", Object_3D_ImageEstimate[i]);
		cv::waitKey(100);
	}

	std::string path = "D://openCV//opencv_file//3DRLD_TEST//RLD//3D RLD//result//";
	//for (int i = 0; i < Nz; i++) {
	//	std::string filename;

	//	cv::Mat img = Object_3D_ImageEstimate[i];
	//	//cv::Mat img_save;
	//	//img.convertTo(img_save, CV_32FC1);
	//	filename = path + "wb_it5"+std::to_string(i + 1) + ".tif";
	//	//cv::imwrite(filename, img_save);
	//	cv::imwrite(filename, img);
	//}
	for (int i = 0; i < Nz; i++) {
		cv::Mat img = Object_3D_ImageEstimate[i]; // CV_64FC1

		//
		cv::Mat save_img;

	

		// 
		cv::Mat normalized;
		double minVal, maxVal;
		cv::minMaxLoc(img, &minVal, &maxVal);

		// 
		const double displayMin = 0.0;
		const double displayMax = 1.0;

		// 
		cv::normalize(img, normalized, displayMin, displayMax, cv::NORM_MINMAX);

		//
		normalized.convertTo(save_img, CV_16UC1, 65535.0);
		cv::imwrite(path + "wb_it2_" + std::to_string(i+1) + ".tif", save_img);
	}



	return 0;
}

