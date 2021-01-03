function real_array = compress_coeff( fft2_coeff )

%fft2 coefficients of real image --> reduced coefficient set

%format of reduced coeff array:
% #c#=one coefficient /// ---c--- multiple coefficients
% RE=real part /// IM=imaginary part
% #RE# ---RE--- #RE# ---IM---
%  |        ---RE---
%  |        ---RE---
%  RE       ---RE---
%  |        ---RE---
%  |        ---RE---
% #RE# ---RE--- #RE# ---IM---
%  |        ---IM---
%  |        ---IM---
%  IM       ---IM---
%  |        ---IM---
%  |        ---IM---


%dimensions
[N,M] = size(fft2_coeff);
N2 = ceil(N/2);
M2 = ceil(M/2);
N_even = (mod(N,2)==0);
M_even = (mod(M,2)==0);

real_array = zeros(N,M);

real_array(1,:) = [real(fft2_coeff(1,1:(M2+M_even))), imag(fft2_coeff(1,(M2+M_even+1):M))];
real_array(2:N2,:) = real(fft2_coeff(2:N2,:));
if N_even == 1
    real_array(N2+N_even,:) = [real(fft2_coeff(N2+N_even,1:(M2+N_even*M_even))), imag(fft2_coeff(N2+N_even,(M2+1+N_even*M_even):M))];
end
real_array((N2+1+N_even):N,:) = imag(fft2_coeff((N2+1+N_even):N,:));

end