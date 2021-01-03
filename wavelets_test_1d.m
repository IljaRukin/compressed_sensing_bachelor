close all
clear
clc

addpath('./wavelets');
addpath('./funktionen');

%Signal
vec=cos(linspace(0,10,128));
%vec=rand(256,1);

fig11=figure(11);set(fig11,'Numbertitle','off','Name','original vector'); clf; plot(vec,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}original vector');

%----- haar
forward_haar = wave_1d_multi(vec,'haar',5);
back_haar = iwave_1d_multi(forward_haar,'haar',5);

%display result
fig1=figure(1);set(fig1,'Numbertitle','off','Name','haar wavelet transform'); clf; plot(forward_haar,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}haar wavelet transform');
fig2=figure(2);set(fig2,'Numbertitle','off','Name','haar reconstructed image'); clf; plot(back_haar,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}haar reconstructed image');

%----- db04
forward_db04 = wave_1d_multi(vec,'db04',5);
back_db04 = iwave_1d_multi(forward_db04,'db04',5);

%display result
fig3=figure(3);set(fig3,'Numbertitle','off','Name','daubechies 4 transform'); clf; plot(forward_db04,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}daubechies 4 transform');
fig4=figure(4);set(fig4,'Numbertitle','off','Name','daubechies 4 reconstructed'); clf; plot(back_db04,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}daubechies 4 reconstructed');

%----- dft and dct
forward_fft = abs( fftshift( fft(vec) * sqrt(numel(vec)) ) );
forward_dct = dct(vec);

%display result
fig5=figure(5);set(fig5,'Numbertitle','off','Name','Absolute value of Fast Fourier Transform'); clf; plot(forward_fft,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}Absolute value of Fast Fourier Transform');
fig6=figure(6);set(fig6,'Numbertitle','off','Name','Discrete Cosine Transform II'); clf; plot(forward_dct,'.');xlabel('t');ylabel('f(t)');title('\fontsize{16}Discrete Cosine Transform II');
