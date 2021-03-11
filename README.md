Methodology for the propagation of a wavefront phase defined in a square-shaped matrix in free space, assuming coherent illumination. At this moment only works for circular pupils.

Developed by Sergio Bonaque-González, PhD.

sergio.bonaque.gonzalez@gmail.com

www.linkedin.com/in/sergiobonaque

July,2019 - Wooptix S.L.

Basically, you define an incoming phase, a propagation distance, and a detector (defined by its pixels number and its length) and you will obtain the intensity image in such a detector.
The software will decide the most appropiate method according to the configuration of the system.

Its an adaptation of the methods described in: 
- Joseph W Goodman - Introduction to Fourier Optics-McGraw-Hill (1996)
- Jason D. Schmidt - Numerical Simulation of Optical Wave Propagation With Examples in MATLAB (2010)
- David Voelz - Computational Fourier Optics (2011)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

It is called in the form of a function of the form:
[x2,Intensity,Uout]=Propagator_ultimate(lambda,z,phase,pupil,L1,L2,N2)

INPUTS:

lambda=wavelength in meters

z=propagation distance in meters

phase=matrix containing the incoming phase defined at the exit pupil 

pupil=matrix containing the exit pupil where the phase is defined

L1=Length of the object space in meters (i.e. length of the matrix containing the wavefront)

L2=Length of the image space in meters (i.e. length of the CCD)

N2=pixels in the detector


OUTPUTS:

x2 = axis values of the Intensity and complex phase results

Intensity = Intensity pattern at the requested z distance

Uout = Complex phase at the requested z distance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ISSUES:
- A comparison with the simpler Fourier transfor method described in Voelz 2011 is also provided beside the final result.

- The first advise "method = xxx" indicates the function used for calculations.

- Next advise indicates the performed interpolation of the incoming phase (bicubic by defect). By definition, twice the resolution is needed as the maximum frequency to be represented. As normally it is not taken into account, it does it automatically. If it has been taken into account when the phase is described, the software should be changed. 

- The last advise indicates the specific method he has used to perform the calculation. This beside the first advise identify the part of the software used for calculations.

- When the direct fourier transform propagation method is used, a geometrical propagation could also be used (not implemented). 

- The boundary conditions to have a good sampling when using the angular
method are very complex to implement. So what I have done is to make sure that the sampling is good enough in the object plane and a warning is displaying saying that if the result is made of several copies or has artifacts, the resolution of the detector must to be increased (normally to the following power of 2).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There is a file called "Example.m" that you can use for testing purposes.

Example of obtained images with the following configuration:

L1=5e-3; % total size of the grid [m]

N1=2^8; %Number of pixels defining the incoming phase

lambda=785e-9; %wavelength [m]

z=0.5; % propagation distance [m]

pupil=mascaraCircular(0.5,N1); %Circular pupil of half the square containing the phase

phase=turbulence_gen(N1).*1e-7; %Generation of a random atmospheric turbulence


%DETECTOR

L2=5e-3; % total size of the grid [m]

N2=2^8;%number of grid points in the detector


![My image1](/imgs/Example_Image.png)   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I am open to include any request or contribution. Do not forget to cite my work if you use it!

