function test2d_1(imgfile,K)
% This example provides a quick comparison of reconstruced images offered
% by the following measurement operators: 
%  1. Scrambled FFT (a dense operator), 
%  2. scrambled 32x32 block Hadamard transform (a highly sparse sampling
%  operator); 
%  3. Block DCT with random sign reversal of input signals (a sparse operator with full streaming capability); 

% Input: 
%   imgfile: input filename;
%   K: No. of Samples; 

% Example: test2d_1('lena256.bmp',25000);
% @Copyright   Kezhi Li 03/12/2009

imgfile='lena256.bmp';
%K=25000;
final=zeros(4,11);

img=imread(imgfile);

% get the row and column dimension of the input image;
[r,c]=size(img);

for jj=1:11
    jj
    K=7000+jj*3000;

result=zeros(4,5);   %change number of loop here
time=zeros(4,1);


for iii=1:5              %%change number of loop here
% Scrambled FFT:
tic;
[rec_sfft,psnr_sfft]=fast_cs2d_frost_1(imgfile, K, 'FFT', 0,32);
result(1,iii)=psnr_sfft;
time(1)=toc+time(1);
% Scrambled Block Walsh-Hadamard Transform;
tic;
[rec_sbhe,psnr_sbhe]=fast_cs2d_frost_1(imgfile, K, 'BWHT', 0,32);
result(2,iii)=psnr_sbhe;
toc;
time(2)=toc+time(2);
% Block DCT with random sign reversal;
% Block size=row number of input image;
tic;
[rec_bdct,psnr_bdct]=fast_cs2d_frost_1(imgfile, K, 'BDCT', 0,32);
result(3,iii)=psnr_bdct;
toc;
time(3)=toc+time(3);
%
tic;
[rec_ostm,psnr_ostm]=fast_cs2d_frost_1(imgfile, K, 'OSTM', 0,32);
result(4,iii)=psnr_ostm;
toc;
time(4)=toc+time(4);
close all;
end
final(1,jj)=mean(result(1,:));
final(2,jj)=mean(result(2,:));
final(3,jj)=mean(result(3,:));
final(4,jj)=mean(result(4,:));

end
% Show the reconstructed image and the PSNR values:
figure(1);
imshow(rec_sfft);
message=sprintf('Scrambled FFT, PSNR=%6.2f dB',psnr_sfft);
title(message);

figure(2);
imshow(rec_sbhe);
message=sprintf('Scrambled 32x32 Block WHT, PSNR=%6.2f dB',psnr_sbhe);
title(message);

figure(3);
imshow(rec_bdct);
message=sprintf('DCT + Permutation, PSNR=%6.2f dB',psnr_bdct);
title(message);

figure(4);
imshow(rec_ostm);
message=sprintf('OSTM + Permutation, PSNR=%6.2f dB',psnr_ostm);
title(message);

figure(5);
hold on;
plot(result(1,:),'*');
plot(result(2,:), 'o');
plot(result(3,:), 'x');
plot(result(4,:), '+');
legend('Scrambled FFT','Scrambled 32x32 Block WHT','DCT+Permutation','OSTM+Permutation');
message=sprintf('time1=%6.2f,time2=%6.2f,time3=%6.2f',time(1),time(2),time(3));
xlabel('Circulation');
ylabel('PSNR');
hold off;
time

figure(6);
hold on;
plot([10000/2^16:3000/2^16:40000/2^16],final(1,:),'*-');
plot([10000/2^16:3000/2^16:40000/2^16],final(2,:), 'o');
plot([10000/2^16:3000/2^16:40000/2^16],final(3,:), 'x');
plot([10000/2^16:3000/2^16:40000/2^16],final(4,:), '+');
legend('Scrambled FFT','WHT','DCT','OSTM');
message=sprintf('time1=%6.2f,time2=%6.2f,time3=%6.2f',time(1),time(2),time(3));
xlabel('Sampling Rate(%)');
ylabel('PSNR');
hold off;