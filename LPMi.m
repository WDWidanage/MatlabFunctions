function [Est] = LPMi(x,y,F,fs,poly_order,method)

% Local polynomial method function when the system has a pole at the 
% origin, pure integrator.
%
% Mandatory Inputs:
%   x - Input time signal size N x 1
%   y - Output time signal size N x 1
%
% Optional inputs
%   F - FFT lines size lf x 1
%   fs - Sampling frequency needed to compute ramp spectrum
%   Poly_order - order of local polynomial default value 2, size 1 x 1
%   method.transient - Set to 1 if transients are to be included in the output sprectrum (default method.transient = 1)
%
% Outputs:
%   Est.G(k) -      Estimated FRF with integrator, size lf x 1
%   Est.T(k) -      Estimated Transients, size lf x 1
%   Est.a -         Estimated ramp coefficient
%   Est.Cv(k) -     Noise variance, size lf x 1
%   Est.CG(k) -     FRF variance, size lf x 1
%   Est.Ct(k) -     Transient variance, size lf x 1
%   Est.Ca -        Ramp coefficient variance, size 1 x 1
%   Est.alpha -     polynomial variables size lf x 2*poly_orderdddd
%   Est.Gtilda(k) - Stable linear dynamic estimate G = Gtilda/s 
%   Est.Cg(k) -     Stable FRF (Gtilda) variance, size lf x 1
%
% W.D. Widanage 17/04/2013 (Middle ground)

try, F;
catch,
    N = length(x);
    F = [2:floor(N/2)];
end
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
end
% vectorise inputs
x = x(:);
y = y(:);
F = F(:);


N = length(x); 

% Compute FFTs
X = fft(x)/sqrt(N);
Y = fft(y)/sqrt(N);
X = X(F);
Y = Y(F);

% Initialise
lf = length(F);
% G=zeros(lf,1);
% T=zeros(lf,1);
% a=zeros(lf,1);
% Cv=zeros(lf,1);
% Cg=zeros(lf,1);
% Ct=zeros(lf,1);
% Ca=zeros(lf,1);
% alpha=zeros(lf,2*poly_order);
S = zeros(lf,1);


% Define s vector and divide input
s = 1i*2*pi*fs*[F-1]/N;
X = X./s;


% Using laplace trasform of finite ramp at discrete frequencies define unit
% ramp spectrum
for kk=1:lf
    if F(kk)==1 %Check if dc is present
        S(kk)=N^2/fs/2/sqrt(N);   %Value of ramp sepctrum at dc
    else
        S(kk)=sqrt(-1)*N^2/(2*pi*(F(kk)-1)*fs^2);   
    end
end

if method.transient==1
    ntheta= 2*poly_order+3; %number of parameters
    R=(ntheta+1)/2; % use 2*R+1 frequencies to estimate the ntheta parameters
else
    ntheta= poly_order+2; %number of parameters
    R=ceil((ntheta+1)/2); % use 2*R+1 frequencies to estimate the ntheta parameters
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
            Kn(rr,:)=[X(Ind(rr)), Xr, 1, r_power, S(Ind(rr))];
        else
            Kn(rr,:)=[X(Ind(rr)), Xr, S(Ind(rr))];
        end
        
    end
    
    %Call numerically stable linear least squares function
    [Theta,results] = Lls(Kn,Out);
    
    if method.transient == 1
        %Extract FRF's and their variances
        Est.Gtilda(kk,1) = Theta(1);                  % LTI branch dynamics
        Est.T(kk,1)=Theta(poly_order+2);            % Transient spectrum
        Theta([1,poly_order+2])=[];
        Est.alpha(kk,:)=Theta(1:poly_order*2);      % polynomial variables
        
        Est.Cg(kk,1)=results.paraVar(1);                         % Variance of G
        Est.Ct(kk,1)=results.paraVar(poly_order+2);              % Transient variance
        
        Est.Cv(kk,1)=results.noiseVar;
        
        Est.a(kk,:)=Theta(poly_order*2+1);          % Ramp gradient
        Est.Ca(kk,1)=results.paraVar(ntheta);                    % Gradient variance
    else
        %Extract FRF's and their variances
        Est.Gtilda(kk,1)=Theta(1);                  % LTI branch dynamics
        Est.T(kk,1)=0;                              % Transient spectrum is zero
        
        Est.alpha(kk,:)=Theta(2:poly_order+1);      % polynomial variables
        
        Est.Cg(kk,1)=results.paraVar(1);                         % Variance of G
        Est.Ct(kk,1)=0;                             % Transient variance is zero
        
        Est.Cv(kk,1)=results.noiseVar;
        
        Est.a(kk,:)=Theta(poly_order+2);            % Ramp gradient
        Est.Ca(kk,1)=results.paraVar(poly_order+2);              % Gradient variance
    end
    
    

end

% Divide G by s to obtian frequency response of dynamics with integrator
Est.G = Est.Gtilda./s;

% Divide variance by |s|^2 to obtian variance estimate of frequency response
Est.CG = Est.Cg./(abs(s).^2);           
    