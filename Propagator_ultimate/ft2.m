function G = ft2(g, delta)
%{
function developed by:
- Jason D. Schmidt - Numerical Simulation of Optical Wave Propagation With Examples in MATLAB (2010)

% function G = ft2(g, delta)
%}
G = fftshift(fft2(fftshift(g))) * delta^2;
end
