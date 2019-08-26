%{
Example for the function Propagator_ultimate. It generates a random
atmospheric turbulence, defines the size and sampling in the object plane,
and defines the size and sampling in the image plane (i.e. a detector).

For more documentation, see Propagator_ultimate function.

Developed by Sergio Bonaque-González, PhD.
sergiob@wooptix.com
July,2019 - Wooptix S.L.
%}
close all
clear all

%IMAGE PUPIL PLANE
L1=5e-3; % total size of the grid [m]
N1=2^8;

lambda=785e-9; %wavelength [m]
z=0.005; % propagation distance [m]
pupil=mascaraCircular(0.5,N1);
phase=turbulence_gen(N1).*1e-7;

%DETECTOR
L2=5e-3; % total size of the grid [m]
N2=2^8;%number of grid points in the detector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x2,Intensity,Uout]=Propagator_ultimate(lambda,z,phase,pupil,L1,L2,N2);