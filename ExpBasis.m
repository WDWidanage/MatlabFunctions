function [y,J] = ExpBasis(theta,x)
% Sum of exponenetial basis functions
%
% Input arguments:
%   theta: Vector of coefficients and exponents, nT x 1
%   x: Independent variable, m x 1
%
%  Output arguments:
%   y: Sum of basis exponential basis functions, m x 1
%   J: Jacobian of function, m x nT
%
% W. D. Widanage 10/08/2016 (Awakening!)

nBF = length(theta)/2;              % Number of exponential basis functions
a = theta(1:nBF);                   % Coefficients
tau = theta(nBF+1:2*nBF);           % Time constants
x = x(:);                           % Vectorise
lx = length(x);                     % Number of data points

% Initialise
y = 0;
J1 = nan(lx,nBF);
J2 = nan(lx,nBF);


for ii = 1:nBF
    y = a(ii)*exp(tau(ii)*x) + y;
    
    J1(:,ii) = exp(tau(ii)*x);
    J2(:,nBF+1) = a(ii)*x.*exp(tau(ii)*x);
    
end

% JAcobian matrix
J = [J1,J2];


end

