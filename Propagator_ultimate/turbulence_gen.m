function phasescreen = turbulence_gen(sz)
%{
function developed by:
- Juan M. Trujillo-Sevilla, PhD. Optical/electrical engineer.
trujillo@wooptix.com
%}

% generate the power spectral density values
cx=(-sz:sz);
mx=(ones(2*sz+1,1)*cx).^2;
mr=sqrt(mx+transpose(mx));
psd=0.023*mr.^(-11/3);
psd(sz+1,sz+1)=0;
% generate the random numbers with Gaussian statistics
randomcoeffs=randn(2*sz+1)+1i*randn(2*sz+1);
% phase screen!
phasescreen=real(fft2(fftshift(sqrt(psd).*randomcoeffs)));
phasescreen = imresize(phasescreen,[sz, sz]);
end