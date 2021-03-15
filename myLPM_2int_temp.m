function [Est]= myLPM_2int_temp(X,Y,F,N,fs,poly_order,method)

% My Local polynomial method function when the systems has a pole at the
% origin, pure integrator.

% Inputs:
%   X - FFT of input signal size lf x 1
%   Y - FFT of output signal size lf x 1
%   F - FFT lines size lf x 1
%   N - Length of time signal in samples, needed to compute ramp sepctrum
%   fs - Sampling frequency needed to compute ramp spectrum
%   Poly_order - order of local polynomial default value 2, size 1 x 1
%   method.transient - Set to 1 if transients are to be included in the output sprectrum (default method.transient = 1)
%
% Outputs:
%   Est.G(k) -  Estimated FRF, size lf x 1
%   Est.T(k) -  Estimated Transients, size lf x 1
%   Est.a -     Estimated ramp gradient
%   Est.b -     Estimated quadratic gradient
%   Est.Cv(k) - Noise variance, size lf x 1
%   Est.Cg(k) - FRF variance, size lf x 1
%   Est.Ct(k) - Transient variance, size lf x 1
%   Est.Ca -    gradient variance, size 1 x 1
%   Est.Cb -    quadratic variance, size 1 x 1
%   Est.alpha - polynomial variables size lf x 2*poly_order
%
% W.D. Widanage 11/03/14 (Meh!)

try, poly_order;
catch,
    poly_order=2;
end
try, fs;
catch,
    fs=1;
end
try, method;
catch,
    method.transient=1;
    t=1
end
% vectorise inputs
X=X(:);
Y=Y(:);
F=F(:);

% Initialise
lf=length(F);
G=zeros(lf,1);
T=zeros(lf,1);
a=zeros(lf,1);
Cv=zeros(lf,1);
Cg=zeros(lf,1);
Ct=zeros(lf,1);
Ca=zeros(lf,1);
alpha=zeros(lf,2*poly_order);
S=zeros(lf,1);
S2 = zeros(lf,1);

%Define unit ramp spectrum
L = N/fs;
for kk=1:lf
    if F(kk)==1 %Check if dc is present
        S(kk)=N^2/fs/2/sqrt(N);   %Value of ramp sepctrum at dc
    else
        harm = F(kk)-1;
        S(kk)= sqrt(-1)*N^2/(2*pi*(F(kk)-1)*fs^2);   % Using laplace trasform of finite ramp at discrete frequencies
        S2(kk) = 1i*L^3*(pi*harm-1i)/(2*pi^2*harm^2);
    end
end
% for kk=1:lf
%     if F(kk)==1 %Check if dc is present
%         S(kk)=N^2/fs^2/2;   %Value of ramp sepctrum at dc
%     else
%         dt = 1;
%         s = 1i*2*pi*(F(kk)-1)/(N*dt);
%         S(kk) = 1/s^2;
%     end
% end

if method.transient==1
    ntheta= 2*poly_order+4; %number of parameters
    R = (ntheta)/2; % use 2*R+1 frequencies to estimate the ntheta parameters
else
    ntheta= poly_order+3; %number of parameters
    R = ceil((ntheta)/2); % use 2*R+1 frequencies to estimate the ntheta parameters
end


Rb=R; % use Rb frequencies before kk
Ra=R; % use Ra frequencies after kk



for kk=1:lf
    %Indicies for frequency range
    if (kk<=Rb)
        Ind=[1:2*R+1];
    end
    if (kk>Rb) && (kk<=lf-Ra)
        Ind=[kk-Rb:kk+Ra];
    end
    if (kk>lf-Ra)
        Ind=[lf-2*R:lf];
    end
    %Output column vector
    Out=Y(Ind);
    
    
    Kn=zeros(2*R+1,ntheta);
    
      
    
    for rr=1:2*R+1
        %Taylor series variable r
        r_Taylor=Ind(rr)-kk;
        for ss=1:poly_order             % Taylor series model order
            r_power(1,ss)=r_Taylor^ss;  % r_power for the transient taylor expansion
            Xr=r_power*X(Ind(rr));      % Xr for the frf taylor expansion
        end
        
        
        if method.transient == 1 % estimate with transients
            Kn(rr,:)=[X(Ind(rr)), Xr, 1, r_power, S(Ind(rr)), S2(Ind(rr))];
        else
            Kn(rr,:)=[X(Ind(rr)), Xr, S(Ind(rr)), S2(Ind(rr))];
        end
        
    end
    
    %Call numerically stable linear least squares function
    [Theta,CT,Cv]=Lls(Kn,Out);
    
    if method.transient == 1
        %Extract FRF's and their variances
        Est.G(kk,1)=Theta(1);                       % LTI branch dynamics
        Est.T(kk,1)=Theta(poly_order+2);            % Transient spectrum
        Theta([1,poly_order+2])=[];
        Est.alpha(kk,:)=Theta(1:poly_order*2);      % polynomial variables
        
        Est.Cg(kk,1)=CT(1);                         % Variance of G
        Est.Ct(kk,1)=CT(poly_order+2);              % Transient variance
        
        Est.Cv(kk,1)=Cv;
        
        Est.a(kk,:)=Theta(poly_order*2+1);          % Ramp gradient
        Est.Ca(kk,1)=CT(poly_order*2+1);            % Gradient variance
        
        Est.b(kk,:)= Theta(poly_order*2+2);          % Quadratic gradient
        Est.Cb(kk,1)= CT(poly_order*2+2);            % Quadratic variance
    else
        %Extract FRF's and their variances
        Est.G(kk,1)=Theta(1);                       % LTI branch dynamics
        Est.T(kk,1)=0;                              % Transient spectrum is zero
        
        Est.alpha(kk,:)=Theta(2:poly_order+1);      % polynomial variables
        
        Est.Cg(kk,1)=CT(1);                         % Variance of G
        Est.Ct(kk,1)=0;                             % Transient variance is zero
        
        Est.Cv(kk,1)=Cv;
        
        Est.a(kk,:)=Theta(poly_order+2);            % Ramp gradient
        Est.Ca(kk,1)=CT(poly_order+2);              % Gradient variance
        
        Est.b(kk,:)= Theta(poly_order+2);          % Quadratic gradient
        Est.Cb(kk,1)= CT(poly_order+2);            % Quadratic variance
    end
    
    
    
end