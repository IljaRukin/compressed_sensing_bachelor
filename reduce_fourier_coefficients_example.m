close all
clear
clc

addpath('./wavelets');
addpath('./funktionen');

%signal
pic = 2*(rand(15,14)-0.5);

%discrete fourier transform
fft2_coeff = fft2(pic);

%reduce fourier coefficients (for real signals)
reduced_dft = compress_coeff(fft2_coeff);

%rebuild original fourier coefficients
original_coeff = decompress_coeff(reduced_dft);

%coefficients error between rebuild and original coefficients
difference = abs(fft2_coeff-original_coeff);

fprintf('error: %i \n',max(max(difference)));

%Display
figure('Numbertitle','off','Name','fft_real');imagesc(real(fft2_coeff));axis image;colormap jet;colorbar;caxis([-10 10]);axis off;
figure('Numbertitle','off','Name','fft_imaginary');imagesc(imag(fft2_coeff));axis image;colormap jet;colorbar;caxis([-10 10]);axis off;
figure('Numbertitle','off','Name','compressed_coefficients');imagesc(real(reduced_dft));axis image;colormap jet;colorbar;caxis([-10 10]);axis off;
figure('Numbertitle','off','Name','decompressed_real');imagesc(real(original_coeff));axis image;colormap jet;colorbar;caxis([-10 10]);axis off;
figure('Numbertitle','off','Name','decompressed_imaginary');imagesc(imag(original_coeff));axis image;colormap jet;colorbar;caxis([-10 10]);axis off;
figure('Numbertitle','off','Name','difference');imagesc(difference);axis image;colormap jet;colorbar;caxis([-10 10]);axis off;
