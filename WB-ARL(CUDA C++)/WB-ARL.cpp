#include "RLD.cuh"
#include "WB_back_projector.h"
#include "FourierTransform.h"
#include <iostream> 
#include <cmath> // 
using namespace std;
//#include "kernel.cuh"
//#include "TEST1.h"
const double PI = 3.141592653589793238462643383279;
const double EPS = 1e-8;

//extern "C" int swap_image(cv::cuda::GpuMat src, cv::cuda::GpuMat dst, int h, int w);//
int Nx = 781;          
int	Ny = 781;          
int center_X = 1824;          
int center_Y = 2432;           
int Nz = 76;
//
//
double n_m = 1;
double lambda = 0.55;
double Mag = 40;
double pixelsize = 6.5 / Mag;
double z_step = 0.5;
double NA_obj = 0.75;
double k = 2 * PI / lambda;
std::vector<cv::Mat> AOTF(Nz);
std::vector<cv::Mat>Image_Ori(Nz);
std::vector<cv::Mat>Object_3D_ImageEstimate(Nz);


void RLDeconvolution(const std::vector<cv::Mat>& Object, std::vector<cv::Mat>& PSF,
	std::vector<cv::Mat>& result, int NI) {

	if (Object.empty() || PSF.empty() || NI <= 0) {
		std::cerr << "RLDDeconvolution: 非法输入（图像为空或迭代次数≤0）" << std::endl;
		return;
	}
	int Nz = Object.size();    
	int rows = Object[0].rows; 
	int cols = Object[0].cols; 
	int type = Object[0].type();


	result.clear();
	result.resize(Nz);
	for (int z = 0; z < Nz; z++) {
		result[z].create(rows, cols, type);
		Object[z].copyTo(result[z]); 
	}


	std::vector<cv::Mat> Hpsf(Nz), HI(Nz), Conv(Nz), DV(Nz), DV_Conv(Nz), HDV(Nz), HDV_Conv(Nz), Estimatetemp(Nz), HConv(Nz);
	std::vector<cv::Mat> ImageEstimate = result; 
	std::vector<cv::Mat> Blurred_3D = Object;    


	FourierTransform::ifftshift3(const_cast<std::vector<cv::Mat>&>(PSF)); 
	FourierTransform::fft3(PSF, Hpsf, 1);
	FourierTransform::fftshift3(Hpsf);


	for (int niter = 0; niter < NI; niter++) {
		std::cout << "niter = " << niter << std::endl;


		FourierTransform::fft3(ImageEstimate, HI);
		FourierTransform::fftshift3(HI);

		MatrixOperation::Mat_Mul_cMat3(Hpsf, HI, HConv);
		FourierTransform::ifftshift3(HConv);
		FourierTransform::ifft3(HConv, Conv, 1);

		MatrixOperation::cMat_Div3(Blurred_3D, Conv, DV);
		FourierTransform::fft3(DV, HDV);
		FourierTransform::fftshift3(HDV);

		MatrixOperation::Mat_Mul_cMat3(Hpsf, HDV, HDV_Conv);
		FourierTransform::ifftshift3(HDV_Conv);
		FourierTransform::ifft3(HDV_Conv, DV_Conv, 1);

		MatrixOperation::Mat_Mul3(ImageEstimate, DV_Conv, Estimatetemp);
		ImageEstimate = Estimatetemp;

		// imshow
		int mid_z = Nz / 2; 
		cv::Mat obj_norm, estimate_norm;


		cv::normalize(Object[mid_z], obj_norm, 0.0, 1.0, cv::NORM_MINMAX, CV_64FC1);

		cv::normalize(ImageEstimate[mid_z], estimate_norm, 0.0, 1.0, cv::NORM_MINMAX, CV_64FC1);


		cv::imshow("Object_3D (Normalized)", obj_norm);
		cv::imshow("RLD - Tra (Normalized)", estimate_norm);
		cv::waitKey(100); 
	}

	result = ImageEstimate;
}

int main()
{
	
	std::vector<cv::Mat>Object(Nz);//
	std::vector<cv::Mat>result(Nz);
	for (int64_t i = 1; i <= Nz; i++)
	{

		/*std::string Image_address = "D://VS_project//RLD//3D RLD//Raw//Rawstack_" + std::to_string(i) + ".tif";*/
		std::string Image_address = "D://openCV//opencv_file//3DRLD_TEST//RLD//3D RLD//Raw//Rawstack_" + std::to_string(i) + ".tif";
		cv::Mat Image_temporary = cv::imread(Image_address, -1); 
		cv::cvtColor(Image_temporary, Image_temporary, cv::COLOR_BGR2GRAY);//
		Image_temporary.convertTo(Image_temporary, CV_64FC1, 1 / 255.0f); //  ->im2double; 8bit->1/255;16bit->1/65535;
		//


		Object[i - 1] = Image_temporary;
		//	cv::imshow("test", Object[i - 1]);
		//cv::waitKey(100);
	}


	double delta_x = 1 / (pixelsize * Nx);
	double delta_y = 1 / (pixelsize * Ny);
	double delta_z = 1 / (z_step * Nz);

	cv::Mat X, Y, Z, fx, fy, fz;
	Matlab2C::linspace_number(-trunc(Nx / 2) * pixelsize, trunc((Nx - 1) / 2) * pixelsize, Nx, X); // linespace X Y Z spatial coordinate
	Matlab2C::linspace_number(-trunc(Ny / 2) * pixelsize, trunc((Ny - 1) / 2) * pixelsize, Ny, Y);
	Matlab2C::linspace_number(-trunc(Nz / 2) * z_step, trunc((Nz - 1) / 2) * z_step, Nz, Z);

	Matlab2C::linspace_number(-trunc(Nx / 2) * delta_x, trunc((Nx - 1) / 2) * delta_x, Nx, fx); // linespace fx fy fz frequency coordinate
	Matlab2C::linspace_number(-trunc(Ny / 2) * delta_y, trunc((Ny - 1) / 2) * delta_y, Ny, fy);
	Matlab2C::linspace_number(-trunc(Nz / 2) * delta_z, trunc((Nz - 1) / 2) * delta_z, Nz, fz);

	cv::Mat Sx, Sy, Fx, Fy;
	MatrixOperation::meshgrid(X, Y, Sx, Sy);
	MatrixOperation::meshgrid(fx, fy, Fx, Fy);

	std::vector<cv::Mat> rhox(Nz), rhoy(Nz), eta(Nz), etaf(Nz);
	cv::Mat MASKeta;
	cv::Mat eta1;

	std::cout << "caculate F" << std::endl;
	for (int zz = 0; zz < Nz; zz++)
	{
		eta[zz] = cv::Mat::ones(cv::Size(Nx, Ny), CV_64FC1);
		eta1 = cv::Mat::ones(cv::Size(Nx, Ny), CV_64FC1);
		eta1 = eta1 * fz.at<double>(zz);
		MASKeta = eta1 == 0;
		MASKeta.convertTo(MASKeta, CV_64FC1, 1 / 255.0f);//double

		eta[zz] = eta1 + (eta1 + EPS).mul(MASKeta);
		etaf[zz] = -eta[zz];
		rhox[zz] = Fx;
		rhoy[zz] = Fy;
	}

	double Max_frequency = NA_obj / lambda; // cutoff frequency
	double coh_para1 = 1; // coherence degree
	double rhos = Max_frequency * coh_para1;
	double rhop = Max_frequency;

	std::vector<cv::Mat> rhoxy(Nz), rho(Nz);
	for (int zz = 0; zz < Nz; zz++)
	{
		cv::sqrt(Fx.mul(Fx) + Fy.mul(Fy), rhoxy[zz]);
	}
	for (int zz = 0; zz < Nz; zz++)
	{
		cv::sqrt(rhoxy[zz].mul(rhoxy[zz]) + eta[zz].mul(eta[zz]), rho[zz]);
	}

	std::vector<cv::Mat> Fp(Nz), Fm(Nz);
	Fp = TransferFunction::Caculate_F(rhoxy, eta, lambda, rho, rhos, rhop);
	Fm = TransferFunction::Caculate_F(rhoxy, etaf, lambda, rho, rhos, rhop);


	for (int zz = 0; zz < Nz; zz++)
	{
		AOTF[zz] = (Fp[zz] + Fm[zz]) * lambda / PI / 4;

	}

	for (int zz = 0; zz < Nz; zz++)
	{
		for (int yy = 0; yy < Ny; yy++)
		{
			for (int xx = 0; xx < Nx; xx++)
			{
				double pv = AOTF[zz].at<double>(yy, xx);
				if (isnan(pv)) {
					AOTF[zz].at<double>(yy, xx) = 0; 
				}
			}
		}
	}
	
	std::vector<cv::Mat> PSF(Nz);
	FourierTransform::ifftshift3(AOTF);
	FourierTransform::ifft3(AOTF, PSF, 1);
	FourierTransform::fftshift3(PSF);

	for (int zz = 0; zz < Nz; zz++)
	{
		MatrixOperation::Nonnegative(PSF[zz]);
	}
	std::cout << "complete caculate F" << std::endl;

	std::vector<cv::Mat> flippedPSF(Nz);

	flipPSF(PSF, flippedPSF);




std::string bp_type = "wiener-butterworth";


std::vector<cv::Mat> PSF_bp = calculatePSF_bp(
	flippedPSF,
	bp_type,
	0.001,   // alpha
	0.01,    // beta
	10,      // n_order
	{ 3.9751, 3.9751, 1.9827 } // iRes
);
	

	WienerButterworth_ARLDeconvolution(Object, PSF, PSF_bp, result ,2);
	NativeRLDeconvolution(Object, AOTF, result,200);
	
	//RLDeconvolution(Object, PSF, result, 200);
	return 0;

}





	//
	
	//

