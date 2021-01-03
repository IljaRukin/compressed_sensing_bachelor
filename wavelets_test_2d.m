close all
clear
clc

addpath('./wavelets');
addpath('./funktionen');

%Bild
pic = struct2array(load('./bilder/dots_1024.mat'));
pic = normalize(pic); %Skalierung auf Werte zwischen -1 und 1

fig11=figure(11);set(fig11,'Numbertitle','off','Name','original picture');imagesc(pic);title('\fontsize{16}original picture');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);

%----- haar

forward_standard = wave_2d_standard(pic,'haar',[1,2]);
%forward=forward.*(abs(forward)>10); %compression
backward_standard = iwave_2d_standard(forward_standard,'haar',[1,2]);
forward_nonstandard = wave_2d_nonstandard(pic,'haar',2);
%forward=forward.*(abs(forward)>10); %compression
backward_nonstandard = iwave_2d_nonstandard(forward_nonstandard,'haar',2);

%display result
fig1=figure(1);set(fig1,'Numbertitle','off','Name','separable haar wavelet transform');imagesc(edit4plot(forward_standard));title('\fontsize{16}separable haar wavelet transform');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig2=figure(2);set(fig2,'Numbertitle','off','Name','separable haar wavelet reconstructed');imagesc((backward_standard));title('\fontsize{16}separable haar wavelet reconstructed');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig3=figure(3);set(fig3,'Numbertitle','off','Name','nested haar wavelet transform');imagesc(edit4plot(forward_nonstandard));title('\fontsize{16}nested haar wavelet transform');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig4=figure(4);set(fig4,'Numbertitle','off','Name','nested haar wavelet reconstructed');imagesc((backward_nonstandard));title('\fontsize{16}nested haar wavelet reconstructed');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);

%----- db04

forward_standard = wave_2d_standard(pic,'db04',[2,2]);
%forward=forward.*(abs(forward)>10); %compression
backward_standard = iwave_2d_standard(forward_standard,'db04',[2,2]);
forward_nonstandard = wave_2d_nonstandard(pic,'db04',2);
%forward=forward.*(abs(forward)>10); %compression
backward_nonstandard = iwave_2d_nonstandard(forward_nonstandard,'db04',2);

%display result
fig5=figure(5);set(fig5,'Numbertitle','off','Name','separable daubechies 4 wavelet transform');imagesc(edit4plot(forward_standard));title('\fontsize{16}separable daubechies 4 wavelet transform');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig6=figure(6);set(fig6,'Numbertitle','off','Name','separable daubechies 4 reconstructed');imagesc((backward_standard));title('\fontsize{16}separable daubechies 4 reconstructed');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig7=figure(7);set(fig7,'Numbertitle','off','Name','nested daubechies 4 wavelet transform');imagesc(edit4plot(forward_nonstandard));title('\fontsize{16}nested daubechies 4 wavelet transform');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig8=figure(8);set(fig8,'Numbertitle','off','Name','nested daubechies 4 reconstructed');imagesc((backward_nonstandard));title('\fontsize{16}nested daubechies 4 reconstructed');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);

%----- fff & dct

forward_fft = fftshift( fft2(pic) / sqrt(numel(pic)) );
forward_dct = dct2(pic);

fig9=figure(9);set(fig9,'Numbertitle','off','Name','Absolute value of Fast Fourier Transform');imagesc(edit4plot(forward_fft));title('\fontsize{16}Absolute value of Fast Fourier Transform');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);
fig10=figure(10);set(fig10,'Numbertitle','off','Name','Discrete Cosine Transform II');imagesc(edit4plot(forward_dct));title('\fontsize{16}Discrete Cosine Transform II');colormap('jet');c = colorbar;c.FontSize = 12;axis off;caxis([-8 2]);


function [output] = edit4plot(input)
    output = log10(abs(input));
end