function [G, T, Cv ,Cg, Ct, alpha] = LPM(X,Y,F,varargin)
%
% Local polynomial method function to estimate a frequency response
% function
% Madotory inputs:
%   X - FFT of input signal size lf x 1
%   Y - FFT of output signal size lf x 1
%   F - FFT lines size lf x 1
%
% Optional arguments:
%   method.order - order of local polynomial default value 2, size 1 x 1
%   method.transient - Set to 1 if transients are to be included in the output sprectrum (default method.transient = 1)
%   method.plot - Set to 1 to plot estimated G and std, default 0
%
% Outputs:
%   G - estimated FRF, size lf x 1
%   T - Estimated Transients, size lf x 1
%   Cv - Noise variance, size lf x 1
%   Cg - FRF variance, size lf x 1
%   Ct - Transient variance, size lf x 1
%   alpha - polynomial variables size lf x 2*poly_order
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 28/06/11-14/11/2020 (Through the never)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments
X = X(:);
Y = Y(:);
F = F(:);

% Create variable names and assign default values after checking the value
addRequired(parObj,'X', @isnumeric);
addRequired(parObj,'Y', @isnumeric);
addRequired(parObj,'F', @isnumeric);


% Optional parameteres
addParameter(parObj,'order',2);
addParameter(parObj,'transient',1);
addParameter(parObj,'plot',0);


% Re-parse parObj
parse(parObj,X,Y,F,varargin{:})

% vectorise inputs
X = parObj.Results.X;
Y = parObj.Results.Y;
F = parObj.Results.F;

% Initialise
lf = length(F);
G = zeros(lf,1);
T = zeros(lf,1);
Cv = zeros(lf,1);
Cg = zeros(lf,1);
Ct = zeros(lf,1);

poly_order = parObj.Results.order;
transient = parObj.Results.transient;
plotLPM = parObj.Results.plot;

if transient == 1
    alpha=zeros(lf,2*poly_order);
else
    alpha=zeros(lf,poly_order);
end

if transient == 1
    ntheta= 2*poly_order+2;
else
    ntheta= poly_order+1;
end

% Use 2*R+1 frequencies to estimate the ntheta parameters
if mod(ntheta,2)==0
    R = ntheta/2;
else
    R = (ntheta + 1)/2;
end

Rb = R; % use Rb frequencies before kk
Ra = R; % use Ra frequencies after kk

for kk = 1:length(F)
    %Slice input and output over frequency range
    if (kk<=Rb)
        Fslice = F(1:2*R+1);
        Xslice = X(1:2*R+1);
        Yslice = Y(1:2*R+1);
    end
    if (kk>Rb) && (kk<=lf-Ra)
        Fslice = F(kk-Rb:kk+Ra);
        Xslice = X(kk-Rb:kk+Ra);
        Yslice = Y(kk-Rb:kk+Ra);
    end
    if (kk>lf-Ra)
        Fslice = F(lf-2*R:lf);
        Xslice = X(lf-2*R:lf);
        Yslice = Y(lf-2*R:lf);
    end
    
    % Regressor matrix
    Kn = zeros(2*R+1,ntheta);
    
    for rr = 1:2*R+1
        %Taylor series variable r
        r_Taylor = Fslice(rr)-F(kk);
        for ss = 1:poly_order             % Taylor series model order
            r_power(1,ss) = r_Taylor^ss;  % r_power for the transient taylor expansion
            Xr = r_power*Xslice(rr);      % Xr for the frf taylor expansion
        end
        if transient == 1 % estimate with transients
            Kn(rr,:)=[Xslice(rr), Xr, 1, r_power];
        else
            Kn(rr,:)=[Xslice(rr), Xr];
        end
    end
    
    %Call numerically stable linear least squares function
    [Theta,results] = Lls(Kn,Yslice);
    
    %Extract FRF's and their variances
    G(kk,1) = Theta(1);                     % LTI branch dynamics
    Cg(kk,1) = results.paraVar(1);          % Variance of G
    Cv(kk,1) = results.noiseVar;            % Noise variance
    
    if transient == 1
        T(kk,1) = Theta(poly_order+2);      % Transient spectrum
        Theta([1,poly_order+2])=[];
        alpha(kk,:) = Theta;                % polynomial variables
        Ct(kk,1) = results.paraVar(poly_order+2);    % Transient variance
    else
        T(kk,1) = nan;                      % Transient spectrum
        Theta(1) = [];
        alpha(kk,:) = Theta;                % polynomial variables
        Ct(kk,1) = nan;                     % Transient variance
    end
end

db = @(x) 20*log(abs(x));

if plotLPM
    figure
    phG = @(G) unwrap(angle(G))*180/pi;     % Function for transfer function phase
    ax1 = subplot(2,1,1);
    semilogx(F,db(G),'. -',F,db(Cg)/2,'. -')
    xlabel ('F [index]'); ylabel('G and stdG [dB]')
    ax2 = subplot(2,1,2);
    semilogx(F,phG(G),'. -')
    xlabel ('F [index]'); ylabel('Phase [deg]')
    linkaxes([ax1,ax2],'x')
end
