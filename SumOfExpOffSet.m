function [y,J] = SumOfExpOffSet(theta,x)
% Sum of exponentials with bias y = U + sum_i^M aexp(bx).
%
% Input arguments:
%   theta: [U, a1,...,aM,b1,...,bM], size M+1 x 1 with M the number of
%          exponentials
%       x: Exponent values, size N x 1
%
% Output
%   y: function value size N x 1
%   J: Jacobian matrix size N x M+1
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 25/12/2016 (Heavy!)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only


p = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
addRequired(p,'theta', @isnumeric);
addRequired(p,'x', @isnumeric);

% Re-parse parObj
parse(p, theta, x)

x = p.Results.x(:);
theta = p.Results.theta(:);

nData = length(x);

nExp = length(theta(2:end))/2; % Number of exponentials
U = theta(1);
a = theta(2:nExp+1);
b = theta(nExp+2:end);

y = zeros(nData,1);

for ii = 1:nExp
    yTmp = a(ii)*exp(b(ii)*x);
    y = y + yTmp;
    J1(:,ii) = exp(b(ii)*x);
    J2(:,ii) = a(ii)*x.*exp(b(ii)*x);
end
y = y + U;
J = [ones(nData,1),J1,J2];


end

