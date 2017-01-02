function [Gsplit, Jsplit] = TFfrf(theta,w,nb,na,varargin)
%
% Evaluates the frequency response and Jacobian of a rational rtransfer
% function. The real and imaginary parts of the functions are stacked and
% returned as outputs. TFfrf is intended for use in transfer optimisation
% routine.
% Mandotory input arguments
%   theta: Numerator and denomentaor coefficients of transfer. Aranaged in
%          decreasing order. theta = [b_nb;b_nb-1;...;b_0;a_na;a_na-1;...;a_1]
%          Coefficient a_0 is fixed to 1. Size  = ntheta x 1
%      w : Vector of normalised angular frequencies rad/sample. Size nw x 1
%     nb : Numeraotr order
%     na : Denomenator order
%
% Optional input arguments. Create a structure variable with the following fields:
%   domain: 's' or 'z' for continous or a discrete transfer fucntion
%       fs: Sampling frequency (Hz) for continous transfer fucnction
%
% Outputs:
%      Gsplit: FRF split into real and imaginary parts
%      Jsplit: Jacobian matrix spilt into real and imaginary parts
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 18/01/2016 (The Hills!!)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only

pObj = inputParser; % Create an input parse object to handle positional and property-value arguments
theta = theta(:);
w = w(:);
nW = length(w);
nPara = length(theta);

% Create required variable names
addRequired(pObj,'theta', @isnumeric);
addRequired(pObj,'w', @isnumeric);
addRequired(pObj,'nb', @isnumeric);
addRequired(pObj,'na', @isnumeric);

% Optional parameteres
addParameter(pObj,'domain','s')
addParameter(pObj,'fs',1)

% Re-parse parObj
parse(pObj,theta,w,nb,na,varargin{:})

bCoeff = theta(1:pObj.Results.nb+1);            % Numerator coefficients
aCoeff = [theta(pObj.Results.nb+2:nPara);1];    % Denomenator coefficients

if pObj.Results.domain =='s' % Continous transfer function
    v = 1i * pObj.Results.w * pObj.Results.fs;
else
    v = exp(1i*pObj.Results.w);
end

bPoly = polyval(bCoeff,v);
aPoly = polyval(aCoeff,v);
J = zeros(nW,nPara);
for kk = 1:nW
    Bkk = bPoly(kk);
    Akk = aPoly(kk);
    dGdb = (1/Akk)*(v(kk).^(nb:-1:0));
    dGda = (-Bkk/Akk^2)*(v(kk).^(na:-1:1));
    J(kk,:) = [dGdb,dGda];
end

G = bPoly./aPoly;

Gsplit = [real(G);imag(G)];
Jsplit = [real(J);imag(J)];

