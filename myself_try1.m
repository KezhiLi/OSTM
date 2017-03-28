imgfile='lena256.bmp';
K=25000;

img=imread(imgfile);

% get the row and column dimension of the input image;
[r,c]=size(img);

% Scrambled FFT:

%[rec_sfft,psnr_sfft]=fast_cs2d_frost_1(imgfile, K, 'FFT', 0);
%function [rec_img,psnr_val]=fast_cs2d_frost_1(im, K, trans_mode, rand_type, blk_size);
im=imgfile;
trans_mode='FFT';
rand_type=0;
blk_size=8;

x = double(imread(im));
x0 = x;

[m n] = size(x);
% Total number of pixels;
N = m*n;
x = x(:);
% Substract the mean of the input signal; 
xmean = mean(x);
x = x-xmean;

% Sparsifying transform: 9-7 Wavelet 
[h0,h1,f0,f1] = filter9_7();
L = floor(log2(m))-3;               % Level of decomposition

% Initialize the output parameters;
rec_img=[];
psnr_val=[];

% Choose an arbitrary random seed; 
% User can change it
Perm_state = 3587642;
rand('state', Perm_state);
  
% Define the random vector for the input samples: 

% other modes: randomize samples use permuation vector
        rand_vect = randperm(N)';
Ki=1;
% Define selected samples
            select_vect = randperm(round(N/2)-1)+1;
            select_vect = select_vect(1:round(Ki/2))';
            % Define Sampling Operator;
            Phi = @(z) fft1d_f(z, select_vect, rand_vect);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) fft1d_t(z, N, select_vect, rand_vect);  

   B = @(z) A_idwt1d(Phi,z,f0,f1,L,m,n);
   Bt = @(z) At_dwt1d(Phi_T,z,h0,h1,L,m,n);
   

    % getting measurements
     y = Phi(x);
     
    % Reconstruction. Using GPSR_BB modules
    tau = norm(y,'fro')/sqrt(Ki)/16;   
    [alp,alp_debias,objective,times,debias_start,mses]= ...
         GPSR_BB(y,B,tau,...
         'AT', Bt,'Debias',1,'Initialization',0,...
         'StopCriterion',1,'ToleranceA',0.001,'ToleranceD',0.00005);
     
   % Transform from the WT domain to the spatial domain  
    alp_debias = reshape(alp_debias,m,n);
    xr = idwt2d(alp_debias,f0,f1,L);
    % add mean to the reconstruction
    xr = xr+xmean;
    psnr_val  = [psnr_val,psnr(x0,xr)];

if length(K)==1,
    rec_img=uint8(xr);
end


rec_sfft=rec_img;
psnr_sfft=psnr_val;



figure(1);
imshow(rec_sfft);
message=sprintf('Scrambled FFT, PSNR=%6.2f dB',psnr_sfft);
title(message);