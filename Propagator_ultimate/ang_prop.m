function [x2, Uout] = ang_prop(config,object,image)
%{
Its an adaptation of the methods described in: 
- Joseph W Goodman - Introduction to Fourier Optics-McGraw-Hill (1996)
- Jason D. Schmidt - Numerical Simulation of Optical Wave Propagation With Examples in MATLAB (2010)
- David Voelz - Computational Fourier Optics (2011)

Adapted by Sergio Bonaque-Gonzalez, PhD. Optical Engineer
sergiob@wooptix.com

July,2019 - Wooptix S.L.

[x2, Uout] = ang_prop(config,object,image)
%}
% SAMPLING CONSTRAINS IN OBJECT PLANE
Nmin= object.L/(2*object.delta) + image.L/(2*image.delta) + (config.lambda*config.z)/(2*object.delta*image.delta);
N_Angular = 2^ceil(log2(Nmin));

Nmin_pre=Nmin+100;
n=1;

while N_Angular>config.resolution_limit && (Nmin_pre-Nmin)>1 %Second condition indicates when further partial propagations does not significally improve the condition
    Nmin_pre=Nmin;
    dz=config.z;
    dz=dz/n;
    Nmin= object.L/(2*object.delta) + image.L/(2*image.delta) + (config.lambda*dz)/(2*object.delta*image.delta);
    N_Angular = 2^ceil(log2(Nmin));
    n=n+1;
end

if N_Angular>config.resolution_limit
    N_Angular=config.resolution_limit;
    fprintf('Reached defined resolution limit=%i. Increase it or check the initial conditions for better sampling.\n',N_Angular);
end

FresnelNumber=object.L^2/(config.z*config.lambda);
if FresnelNumber>1500 && object.delta~=image.delta %Value obtained experimentally
    if length(object.phase)<N_Angular
        fprintf('Incoming phase has been resized to %ix%i because sampling constrains in object plane.\n',N_Angular,N_Angular);
        
        [first, last]= get_fist_last_non_zero_index(object.pupil);
        pupil=mascaraCircular((last-first)/length(object.pupil), N_Angular);
        phase=imresize(object.phase, [N_Angular N_Angular]);
        object.delta=object.L/N_Angular;
    else
        pupil=object.pupil;
        phase=object.phase;
    end
    Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
    [~, Uout]= propTF(Uin,object.L,config.lambda,config.z);
    fprintf('F-NUMBER EXTREMELY HIGH. DIRECT FOURIER TRANSFORM PROPAGATION. SIMPLE INTERPOLATION TO FIT THE DETECTOR\n'); %valid for very short propagations
    x2=-object.L/2:object.delta:object.L/2-object.delta;
    if object.L>image.L
        for i=1:length(x2)
            if x2(i)>=-image.L/2
                first=i;
                break;
            end
        end
        last=length(x2);
        for i=1:length(x2)
            if x2(i)>=image.L/2
                last=i;
                break;
            end
        end
        x2=x2(first:last);
        Uout=Uout(first:last,first:last);
        x2=imresize(x2,[1 image.N]);
        Uout=imresize(Uout,[image.N image.N]);
    elseif object.L<image.L
        n=1;
        while object.L+n*object.delta<image.L
            n=n+1;
        end
        if rem(n,2)~=0
            n=n+1;
        end
        Uout = padarray(Uout,[n/2 n/2],0,'both');
        Uout=imresize(Uout,[image.N image.N]);
        x2=-image.L/2:image.delta:image.L/2-image.delta;
    end
    
    
    
else
    if length(object.phase)<N_Angular
        fprintf('Incoming phase has been resized to %ix%i because sampling constrains in object plane.\n',N_Angular,N_Angular);
        [first, last]= get_fist_last_non_zero_index(object.pupil);
        pupil=mascaraCircular((last-first)/length(object.pupil), N_Angular);
        phase=imresize(object.phase, [N_Angular N_Angular]);
        object.delta=object.L/N_Angular;
    else
        pupil=object.pupil;
        phase=object.phase;
    end
    if n==1
        fprintf('--Angular spectrum method--\nThe analysis of proper resolution in the detector for angular propagation is complex.\nWARNING:If multiple copys or artifacts are clearly visible, increase virtually the resolution of the detector.\n'); %valid for short propagations
        Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
        
        [x2, Uout] = ang_spec_prop(Uin, config.lambda, object.delta, image.delta,config.z);
    else
        fprintf('--Angular spectrum method with %i partial propagations--\nThe analysis of proper resolution in the detector for angular propagation is complex.\nWARNING: If multiple copys or artifacts are clearly visible, increase virtually the resolution of the detector.\n',n); %valid for short propagations
        z = (1:n) * config.z / n;
        Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
        [x2, Uout] = ang_spec_multi_prop_vac (Uin, config.k, object.delta, image.delta, z);
        
    end
    
    
    if image.N<length(Uout)
        m=floor((length(Uout)-image.N)/2);
        Uout(end-m+1:end,:)=[];
        Uout(:,end-m+1:end)=[];
        Uout(1:m,:)=[];
        Uout(:,1:m)=[];
        x2(end-m+1:end,:)=[];
        x2(:,end-m+1:end)=[];
        x2(1:m,:)=[];
        x2(:,1:m)=[];
        if rem(image.N,2)~=0
            Uout(:,end)=[];
            Uout(end,:)=[];
            x2(:,end)=[];
            x2(end,:)=[];
        end
    end
end
