function pupil = mascaraCircular(radius,resolution)
%{
Creates a pupil with a defined radius (in %) inside a square matrix of
given resolution

by Sergio Bonaque-Gonzalez, PhD. Optical Engineer
sergiob@wooptix.com
July,2019 - Wooptix S.L.
%}
xp=linspace(-1,1,resolution);
[X,Y]=meshgrid (xp,xp); 
[rho]=sqrt(X.^2+Y.^2); 

pupil=ones(size(rho));
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > radius
            pupil(i,j)=0;
        end;
    end;
end;

