close all
clear
clc

addpath('./wavelets');
addpath('./funktionen');
addpath('./solver/implemented');

%load image
pic = struct2array(load('./bilder/dots_256.mat'));
pic = normalization(pic);

%display image
figure(1);imagesc(pic);colormap(jet);colorbar();

%generate samples
dim = size(pic);
mask = mask(dim,0.1);
indices = find(mask);
pic_samples = pic(indices);

%display samples
punkte = zeros(dim);
punkte(indices) = pic_samples;
figure(4);imagesc(punkte);colormap(jet);colorbar;

%select transformations (basis)
trafo = 'cdft';
if strcmp(trafo,'dft2')
    AA = @(x) fft2(x)/sqrt(prod(dim));
    AT = @(x) ifft2(x)*sqrt(prod(dim));
elseif strcmp(trafo,'dct2')
    AA = @(x) dct2(x);
    AT = @(x) idct2(x);
elseif strcmp(trafo,'cdft')
    AA = @(x) compress_coeff(fft2(x))/sqrt(prod(dim));
    AT = @(x) ifft2(decompress_coeff(x))*sqrt(prod(dim));
elseif strcmp(trafo,'haar')
    divisibility = 0;divider=2;
    while rem(dim,divider)==0
        divisibility = divisibility +1; divider = divider*2;
    end
    %%%
    AA = @(x) wave_2d_nonstandard(x,'haar',divisibility);
    AT = @(x) iwave_2d_nonstandard(x,'haar',divisibility);
elseif strcmp(trafo,'db04')
    divisibility = 0;divider=2;
    while rem(dim,divider)==0
        divisibility = divisibility +1; divider = divider*2;
    end
    %%%
    AA = @(x) wave_2d_nonstandard(x,'db04',divisibility);
    AT = @(x) iwave_2d_nonstandard(x,'db04',divisibility);
else
    fprintf('wrong transformation basis ! \n chose one of the following: spgl1 nesta fpcas');
end

%OMP parameters
sparsity = 200; %maximal number of coefficients. When reached, algorithm stops execution.
start_lsqlin_OptimalityTolerance = 1e-1; %accuracy for solving each coefficient value
final_lsqlin_OptimalityTolerance = 1e-3; %accuracy for solving final coefficients values
accuracy = 0.1; %fit error for samples. When reached, algorithm stops execution.

%%%%% normalize image: mean=0 & variance=1
mean = sum(sum(sum(pic_samples)))/numel(pic_samples);
pic_samples = pic_samples - mean;
variance = sum(sum(sum(pic_samples.^2)))/numel(pic_samples);
pic_samples = pic_samples/sqrt(variance);

%start solver
tic
x = solveOMP(pic_samples,dim,indices,sparsity,start_lsqlin_OptimalityTolerance,accuracy,final_lsqlin_OptimalityTolerance,AA,AT);
comp_time = toc;
result = AT(x);
result = real(result);

%reverse normalization
result = result*sqrt(variance) + mean;
fit_error = sum(sum(sum((pic_samples-result(indices)).^2)))/sum(sum(sum(pic_samples.^2)));

%display result
fprintf('computation time: %f \n',comp_time);
fprintf('image error: %i \n',sum(sum(sum((pic-result).^2)))/sum(sum(sum(pic.^2))));
fprintf('fitting error: %i \n',fit_error);
fprintf('l1 norm image: %f \n',sum(sum(sum(abs(pic)))));
fprintf('l1 norm reconstructed: %f \n',sum(sum(sum(abs(result)))));
figure(2);imagesc(result);title('reconstructed');colormap(jet);colorbar;
figure(3);imagesc(result-pic);title('error');colormap(jet);colorbar;
