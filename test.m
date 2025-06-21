%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for Experiment for WB-ARL.
% The Code is created based on the method described in the following paper 
%   [1]Du H, Zhang R, Zuo C, et al.
%      "Wiener–Butterworth Accelerated Richardson–Lucy Deconvolution for Rapid 3D Fluorescence Microscopy"
%   [2]Guo, Min, et al. "Rapid image deconvolution and multiview fusion for optical microscopy." Nature biotechnology 38.11 (2020): 1337-1346. 
%      The Code is modified and extended from Guo's code
%    Contact: HeHeng Du (hehengdu@njust.edu.cn)
%    Date  : 2025/6/20

clc
clear 
close all
tic
set(0, 'DefaultFigureWindowStyle', 'docked');%

 

path = './Raw';
A = im2double(rgb2gray(imread([path,'/Rawstack_35','.tif'])));
figure,imshow(A,[])

Zn=0;
for num = 1:1:76
    Zn = Zn+1;
    Z_array(:,:,Zn)= im2double(rgb2gray(imread([path,'/Rawstack_',int2str(num),'.tif'])));

    imshow(Z_array(:,:,Zn),[])
end
[Nx,Ny,Nz]=size(Z_array);


% PSF  读取
input_filename = 'PSF_stack.tif';
info = imfinfo(input_filename);
num_layers = numel(info);

[rows, cols] = size(imread(input_filename, 1));
PSFa = zeros(rows, cols, num_layers, 'uint16');

for z = 1:num_layers
    PSFa(:,:,z) = imread(input_filename, z);
end

load('PSF.mat', 'min_val', 'max_val');
PSFa = double(PSFa) / 65535 * (max_val - min_val) + min_val;

figure,imshow(PSFa(:,:,Nz/2),[]);


%   Generate back_projector
bp_type = 'wiener-butterworth';  
alpha = 0.001;
beta = 0.001; 
n = 10;
resFlag = 1;
iRes = [4.5467,4.5467,2.1526];
verboseFlag = 0;


[PSF_bp, ~] = WB_back_projector(PSFa, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag);


%   WB_back_projector RL Deconvolution Deconvolution 
smallValue = 0;
ImageEstimate_wb = Z_array;
niter=1;
figure,
for i = 1:niter
    disp(i)

    Hpsf = fftshift(fftn(ifftshift(PSFa)));
    HI = (fftshift(fftn(ImageEstimate_wb)));
    Conv = ifftn(ifftshift(Hpsf.*HI));

    DV = Z_array./Conv;         
    HPSFpFlip = fftshift(fftn(ifftshift(PSF_bp)));% PSF_bp:  back_projector
    HDV = fftshift(fftn(DV));
    DV_Conv = ifftn(ifftshift(HDV.*HPSFpFlip));

    ImageEstimate_wb = DV_Conv.*ImageEstimate_wb;
    ImageEstimate_wb = max(ImageEstimate_wb,smallValue);% non-negative constraint
   
    subplot(1,2,1),imshow(Z_array(:,:,Nz/2),[]);title('Blurred')
    subplot(1,2,2),imshow(ImageEstimate_wb(:,:,Nz/2),[]);title('ImageEstimate')
   
end


%   Traditional RL Deconvolution 
ImageEstimate_rld = Z_array;
niter = 20;
figure,
PSF_bp_tra=flipPSF(PSFa);
for i = 1:niter
    disp(i)


    Hpsf = fftshift(fftn(ifftshift(PSFa)));
    HI = (fftshift(fftn(ImageEstimate_rld)));
    Conv = ifftn(ifftshift(Hpsf.*HI));
    DV = Z_array./Conv;     
    HPSFpFlip = fftshift(fftn(ifftshift(PSF_bp_tra)));
    HDV = fftshift(fftn(DV));
    DV_Conv = ifftn(ifftshift(HDV.*HPSFpFlip));


    ImageEstimate_rld = DV_Conv.*ImageEstimate_rld;
    ImageEstimate_rld = max(ImageEstimate_rld,smallValue);     
    
    subplot(1,2,1),imshow(Z_array(:,:,Nz/2),[]);title('Blurred')
    subplot(1,2,2),imshow(ImageEstimate_rld(:,:,Nz/2),[]);title('ImageEstimate')     
end











