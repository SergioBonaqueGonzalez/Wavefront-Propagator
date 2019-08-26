function[fx,u2]=propTF(u1,L,lambda,z)
%{
function developed by:
- David Voelz - Computational Fourier Optics (2011)

% propagation - transfer function approach 
% assumes same x and y side lengths and 
% uniform sampling. This propagator function takes the source field u1 and produces the observation field u2 where the source and observation side lengths and sample coordinates are identical.  
% u1 - source plane field 
% L - source and observation plane side length 
% lambda - wavelength 
% z - propagation distance 
% u2 - observation plane field 
%}

[M,~]=size(u1); %get input field array size 
dx=L/M; %sample interval 

fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords 
[FX,FY]=meshgrid(fx,fx); 
H=exp(-1i*pi*lambda*z*(FX.^2+FY.^2)); %trans func 
H=fftshift(H); %shift trans func 
U1=fft2(fftshift(u1)); %shift, fft src field 
U2=H.*U1; %multiply 
u2=ifftshift(ifft2(U2)); %inv fft, center obs field 
end