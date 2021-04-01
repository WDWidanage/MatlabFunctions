function [mat,info] = ChebLagrangeDerMat(varargin)
%
% Computes the first and second order derivative matrices based on
% 1D Lagrange interpolation functions and Chebyshev collocation points of
% the first kind
%
% Inputs (optional):
%   N: number of collocation points. Double size, 1 x 1. Default 6;
%   int: Domain interval, Double size 1 x 2. Default [-1,1]. Use more
%        collacation points if domain is increased
%
% Outputs:
%   mat: Structure variable with following fields
%        - D: First order derivative matrix. Double size N x N
%        - D2: Second order derivative matrix. Double size N x N
%
%   info: Structure variable with following fields
%            - sn: Collocation points on standard interval [-1,1], size N x 1
%            - xn: Transformed collocation points on to required input interval, size, N x 1
%            - N: Number of collocation points, size 1 x 1
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 15-03-2021 (The Lost Voices)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments

% Optional parameters
addOptional(parObj,'N',6,@isnumeric);
addOptional(parObj,'int',[-1,1],@isnumeric);

% Re-parse parObj
parse(parObj,varargin{:})

N = parObj.Results.N;
int = parObj.Results.int;

P = N - 1;                          % Chebyshev polynomial order
sn = -cos([0:P]'*pi/P);              % Collocation points on standard interval 
xn = int(1) + diff(int)*(1+sn)/2;   % Collocation points on input interval
info.sn = sn;
info.xn = xn;
info.N = N;

% Define derivative matrices
% El-Baghdady et al. Math 3.1 (2016): 1-8.
for ii = 1:N
    for jj = 1:N
        if (ii == jj) && (ii == 1)
            D(1,1) = -(2*P^2+1)/6;
        elseif (ii == jj) && (ii == N)
            D(N,N) = (2*P^2+1)/6;
        elseif (ii == jj) && (ii <= N-1)
            D(ii,jj) = - sn(ii)/(2*(1-sn(ii)^2));
        else
            if (ii == 1) || (ii == N), ci = 2; else, ci = 1; end
            if (jj == 1) || (jj == N), cj = 2; else, cj = 1; end
            D(ii,jj) = ci/cj*(-1)^(ii+jj)/(sn(ii)-sn(jj));
        end
    end
end
mat.D = D*(2/diff(int));     % First order derivative matrix translated to desired domain coordinates. Double size N x N
mat.D2 = D*D;                % Second order derivative matrix tranlated to desired domain coordinates. Double size N x N

end

