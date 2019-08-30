function [x2,Uout]=Fresnel_1step_Propagation(config,object,image)
%{
Its an adaptation of the methods described in: 
- Joseph W Goodman - Introduction to Fourier Optics-McGraw-Hill (1996)
- Jason D. Schmidt - Numerical Simulation of Optical Wave Propagation With Examples in MATLAB (2010)
- David Voelz - Computational Fourier Optics (2011)

This first method evaluates the Fresnel diffraction integral once as a
single FT, which is the most straightforward. It has not control over spacing in the final grid without changing the geometry.
This method is desirable because of its computational efficiency.

Adapted by Sergio Bonaque-Gonzalez, PhD. Optical Engineer
sergiob@wooptix.com
July,2019 - Wooptix S.L.

%}

% SAMPLING CONSTRAINS IN OBJECT PLANE
%The key to achieving an accurate result is to sample the quadratic phase factor inside the FT at a high enough rate to satisfy the Nyquist criterion. If it is not sampled finely enough, the intended high-frequency content would show up in the lower frequencies
% We want to find the maximum local spatial frequency of the quadratic
% phase factor inside the integral and sample at least twice this rate. To
% find it we calculate the maximum slope in the system (from one corner of
% the object plane to the opposite corner of the image plane.
Nmin_Fresnel = object.L * config.lambda*config.z / (object.delta * (config.lambda*config.z - image.L*object.delta));
N_Fresnel = 2^ceil(log2(Nmin_Fresnel));

dz=config.z;
n=1;

Nmin_Fresnel_pre=Nmin_Fresnel+100;
while N_Fresnel>config.resolution_limit && (Nmin_Fresnel_pre-Nmin_Fresnel)>1 %Second condition indicates when further partial propagations does not significally improve the condition
    Nmin_Fresnel_pre=Nmin_Fresnel;
    dz=config.z;
    n=n+1;
    dz=dz/n;
    Nmin_Fresnel = object.L * config.lambda*dz / (object.delta * (config.lambda*dz - image.L*object.delta));
    N_Fresnel = 2^ceil(log2(Nmin_Fresnel));
end

if N_Fresnel>config.resolution_limit
    N_Fresnel=config.resolution_limit;
    fprintf('Reached defined resolution limit=%i. Increase it or check the initial conditions for better sampling.\n',N_Fresnel);
end

%Second criteria, from Voelz
dx=object.L/object.N; %sample interval
crit = abs(config.lambda*(config.z/n)/object.L);
while dx < crit
    n=n+1;
    crit = abs(config.lambda*(config.z/n)/object.L);
end


FresnelNumber=object.L^2/(config.z*config.lambda);
if FresnelNumber>1500 %Value obtained experimentally
    if length(object.phase)<N_Fresnel
        fprintf('Incoming phase has been resized to %ix%i because sampling constrains in object plane.\n',N_Fresnel,N_Fresnel);
        
        [first, last]= get_fist_last_non_zero_index(object.pupil);
        pupil=mascaraCircular((last-first)/length(object.pupil), N_Fresnel);
        phase=imresize(object.phase, [N_Fresnel N_Fresnel]);
        object.delta=object.L/N_Fresnel;
%         m=(N_Fresnel-object.N)/2;
    else
        pupil=object.pupil;
        phase=object.phase;
    end
    Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
    [~, Uout]= propTF(Uin, config.lambda, object.delta, config.z);
    fprintf('--Direct transform propagation. Simple interpolation to fit the detector--\n'); %valid for very short propagations
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
    if dz>=object.L*object.delta/config.lambda
        if n==1
            if length(object.phase)<N_Fresnel
                fprintf('Incoming phase has been resized to %ix%i because sampling constrains in object plane.\n',N_Fresnel,N_Fresnel);
                %             pupil=imresize(object.pupil,[N_Fresnel N_Fresnel],'nearest');
                [first, last]= get_fist_last_non_zero_index(object.pupil);
                pupil=mascaraCircular((last-first)/length(object.pupil), N_Fresnel);
                phase=imresize(object.phase, [N_Fresnel N_Fresnel]);
                object.delta=object.L/N_Fresnel;
                m=floor((N_Fresnel-object.N)/2);
            else
                pupil=object.pupil;
                phase=object.phase;
            end
            
            if object.delta==image.delta
                Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
                fprintf('--Fresnel one-step propagation--\n'); %valid for long propagations
                [x2, Uout]= one_step_prop(Uin, config.lambda, object.delta, config.z);
            else
                Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
                fprintf('--Fresnel two-step propagation--\n'); %valid for long propagations
                [x2, Uout]= two_step_prop(Uin, config.lambda, object.delta, image.delta, dz);
                
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
                    
                else
                    m=length(Uout)/image.N;
                    x2=x2*m;
                end
                
            end
        else
            while dz>=object.L*object.delta/config.lambda %make sure fits with the conditions for angular spectrum propagation
                dz=config.z;
                n=n+1;
                dz=dz/n;
            end
            %make sure fits with sampling constrains
            %Previous methodology can produce very large arrays that can be
            %costly computationally, so the partial propagation method is implemented
            %here in order to maintain the array sizes below or equal to the value
            %specified by config.resolution_limit. So, its sampling constrains are also
            %calculated
            Nmin= object.L/(2*object.delta) + image.L/(2*image.delta) + (config.lambda*dz)/(2*object.delta*image.delta);
            N_Angular = 2^ceil(log2(Nmin));
            Nmin_pre=Nmin+100;
            while N_Angular>config.resolution_limit && (Nmin_pre-Nmin)>1
                Nmin_pre=Nmin;
                dz=config.z;
                dz=dz/n;
                Nmin= object.L/(2*object.delta) + image.L/(2*image.delta) + (config.lambda*dz)/(2*object.delta*image.delta);
                N_Angular = 2^ceil(log2(Nmin));
                n=n+1;
            end
            
            if length(object.phase)<N_Angular
                fprintf('Incoming phase has been resized to %ix%i because sampling constrains in object plane.\n',N_Angular,N_Angular);
                %             pupil=imresize(object.pupil,[N_Angular N_Angular],'nearest');
                [first, last]= get_fist_last_non_zero_index(object.pupil);
                pupil=mascaraCircular((last-first)/length(object.pupil), N_Angular);
                phase=imresize(object.phase, [N_Angular N_Angular]);
                object.delta=object.L/N_Angular;
            else
                pupil=object.pupil;
                phase=object.phase;
            end
            
            
            fprintf('--Angular spectrum method with %i partial propagations--\nThe analysis of proper resolution in the detector for angular propagation is complex.\nWARNING:If multiple copys or artifacts are clearly visible, increase virtually the resolution of the detector.\n',n); %valid for short propagations
            z = (1:n) * config.z / n;
            Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
            [x2, Uout] = ang_spec_multi_prop_vac (Uin, config.k, object.delta, image.delta, z);
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
                
            else
                m=length(Uout)/image.N;
                x2=x2*m;
            end
            
        end
    else
        Nmin= object.L/(2*object.delta) + image.L/(2*image.delta) + (config.lambda*dz)/(2*object.delta*image.delta);
        N_Angular = 2^ceil(log2(Nmin));
        Nmin_pre=Nmin+100;
        while N_Angular>config.resolution_limit && (Nmin_pre-Nmin)>1
            Nmin_pre=Nmin;
            dz=config.z;
            dz=dz/n;
            Nmin= object.L/(2*object.delta) + image.L/(2*image.delta) + (config.lambda*dz)/(2*object.delta*image.delta);
            N_Angular = 2^ceil(log2(Nmin));
            n=n+1;
        end
        
        if length(object.phase)<N_Angular
            fprintf('Incoming phase has been resized to %ix%i because sampling constrains in object plane.\n',N_Angular,N_Angular);
            pupil=imresize(object.pupil,[N_Angular N_Angular],'nearest');
            phase=imresize(object.phase, [N_Angular N_Angular]);
            object.delta=object.L/N_Angular;
        else
            pupil=object.pupil;
            phase=object.phase;
        end
        Uin = pupil.*exp(config.k*1i.*phase); %complex phase screen
        fprintf('--Angular spectrum method with %i partial propagations--\nThe analysis of proper resolution in the detector for angular propagation is complex.\nWARNING:If multiple copys or artifacts are clearly visible, increase virtually the resolution of the detector.\n',n); %valid for short propagations
        z = (1:n) * config.z / n;
        [x2, Uout] = ang_spec_multi_prop_vac (Uin, config.k, object.delta, image.delta, z);
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
            
        else
            m=length(Uout)/image.N;
            x2=x2*m;
        end
    end
end


