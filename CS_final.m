close all
clear
clc

addpath('./wavelets');
addpath('./funktionen');

%load image
pic = struct2array(load('./bilder/dots_1024.mat'));
pic = normalization(pic);

%display image
figure(1);imagesc(pic);colormap(jet),colorbar;

%generate samples
dim = size(pic);
mask = mask(dim,0.1);
indices = find(mask);
pic_samples = pic(indices);

%display samples
samples = zeros(dim);
samples(indices) = pic_samples;
figure(4);imagesc(samples);title('samples');colormap(jet),colorbar;

%set solver and basis
%'cdft','dft2','dct2','haar','db04'
trafo = 'cdft';
%'nesta','spgl1','fpcas'
solver = 'nesta';

%set parameter for noise (reduction)
sigma = 1e-5;

%start solver
[result,fit_error,comp_time] = reconstruct(pic_samples,indices,dim,trafo,solver,sigma);

%display result
fprintf('computation time: %f \n',comp_time);
fprintf('image error: %i \n',sum(sum(sum((pic-result).^2)))/sum(sum(sum(pic.^2))));
fprintf('fitting error: %i \n',fit_error);
fprintf('l1 norm image: %f \n',sum(sum(sum(abs(pic)))));
fprintf('l1 norm reconstructed: %f \n',sum(sum(sum(abs(result)))));
figure(2);imagesc(result);title('reconstructed');colormap(jet);colorbar;
figure(3);imagesc(result-pic);title('error');colormap(jet);colorbar;
