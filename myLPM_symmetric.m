function [G, T, Cv ,Cg, Ct, alpha]= myLPM(X,Y,F,poly_order)

% My Local polynomial method function
% Inputs:
%   X - FFT of input signal size lf x 1
%   Y - FFT of output signal size lf x 1
%   F - FFT lines size lf x 1
%   Poly_order - order of local polynomial default value 2, size 1 x 1
%
% Outputs:
%   G - estimated FRF, size lf x 1
%   T - Estimated Transients, size lf x 1
%   Cv - Noise variance, size lf x 1
%   Cg - FRF variance, size lf x 1
%   Ct - Transient variance, size lf x 1
%   alpha - polynomial variables size lf x 2*poly_order
%
% W.D. Widanage 28/06/11 (Connected?)

try
    if isempty(poly_order);
        poly_order=2;
    end
catch
    poly_order=2;
end

% vectorise inputs
X=X(:);
Y=Y(:);
F=F(:);

% Initialise
lf=length(F);
G=zeros(lf,1);
T=zeros(lf,1);
Cv=zeros(lf,1);
Cg=zeros(lf,1);
Ct=zeros(lf,1);
alpha=zeros(lf,2*poly_order);



ntheta= 2*poly_order+2;
R=ntheta/2; % use R frequencies on either side of kk to estimate the ntheta parameters

for kk=1:lf
    %Indicies for frequency range
    if (kk<=R)
         Ind=[1:2*R+1];
    end
    if (kk>R) & (kk<=lf-R)
        Ind=[kk-R:kk+R];
    end
    if (kk>lf-R)
        Ind=[lf-2*R:lf];
    end
      
    %Output column vector
    Out=Y(Ind);
    
    % Regressor matrix
    Kn=zeros(2*R+1,ntheta);
    for rr=1:2*R+1
        %Taylor series variable r
        r_Taylor=Ind(rr)-kk;
        for ss=1:poly_order             % Taylor series model order           
            r_power(1,ss)=r_Taylor^ss;  % r_power for the transient taylor expansion
            Xr=r_power*X(Ind(rr));      % Xr for the frf taylor expansion
        end
        Kn(rr,:)=[X(Ind(rr)), Xr, 1, r_power];
    end
    
    %Scale regressor for numerical stabilty and perform least squares via
    %SVD; Out=Kn*Theta
    Scale=sqrt(sum(abs(Kn).^2));
    idxZeros = find(Scale<1E-14);
    Scale(idxZeros) = 1;
    Kn_scaled=Kn./repmat(Scale,2*R+1,1);
    %SVD of Kn
    [Un, Sn, Vn] = svd(Kn_scaled,0);
    ss = diag(Sn);
    idxZeros = find(ss < 1E-14);
    ss(idxZeros) = inf;
    ss = diag(1./ss);
    Theta=Vn*ss*Un'*Out;
    Theta=Theta./Scale';    % Theta = [G g1,...,gpoly_order, T, t1,...,tpoly_order]
   
    %Projection matrix P= I-Kn(Kn'Kn)^-1Kn'; use svd of Kn
    Pn=eye(2*R+1)-Un*Un';
    q=real(trace(Pn));
    %Residual
    Rs=Out-Kn*Theta;
    
    %Noise variance estimate;  see pg241 Beck and Arnold
    Cv(kk,1)=(Rs'*Rs)/q;
    
    %Covariance of Theta CT=diag((Kn'Kn)^-1Cv); use svd of Kn see pg239 of
    %Beck and Arnold
    CT=diag(inv(Kn'*Kn)*Cv(kk,1)); % extract variances through diagonal
    
    %Extract FRF's and their variances
    G(kk,1)=Theta(1);                       % LTI branch dynamics
    T(kk,1)=Theta(poly_order+2);            % Transient spectrum
    Theta([1,poly_order+2])=[];  
    alpha(kk,:)=Theta;                      % polynomial variables   
    
    Cg(kk,1)=CT(1);            % Variance of G
    Ct(kk,1)=CT(poly_order+2); % Transient variance
end