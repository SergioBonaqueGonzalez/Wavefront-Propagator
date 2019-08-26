function g = ift2(G, delta_f)
%{
function developed by:
- Jason D. Schmidt - Numerical Simulation of Optical Wave Propagation With Examples in MATLAB (2010)

function g = ift2(G, delta_f)
July,2019 - Wooptix S.L.
%}

N = size(G, 1);
g = ifftshift(ifft2(ifftshift(G))) * (N * delta_f)^2;
end