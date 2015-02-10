function [theta,Ct,Cv]=Lls(K,Z,ny)
% 
% Computes a numerically stable (weighted) linear leaast squares estimate
% and returns the optimal estimate, its variance and an estimate of the
% noise variance
%
% Inputs (mandatory): 
%   K: Regresor matrix, size n x m
%   Z: Output vector, size n x 1
% Inputs (optional): 
%   ny: Output noise vector, std of output noise for weighted least
%   squares, size n x 1. default ny = ones(n,1);
%
% Outputs:
%   theta: Optimum parameter estimate, size m x 1
%   Ct: Covariance estimate of theta, size m x 1
%   Cv: Estimate of noise variance, size 1 x 1
%   
%   Note: Ct and Cv are valid only when the output noise variance is assumed constant
%         and should be disregarded if weighted Lls is performed.
%
% W.D. Widanage 11-02-2012 (Steps)
%

%Check if weight matrix is provided
try, ny;
    %Create weight matrix
    W=diag(1./ny);
catch, 
    W=diag(ones(length(Z),1));
end

%Multiply by weights
Z=W*Z;
K=W*K;

%Normalise K with the 1/norm of each column
Knorm=sqrt(sum(abs(K).^2));
idxZeros = find(Knorm<1E-14);
Knorm(idxZeros) = 1;
N=diag(1./Knorm);
Kn=K*N;

%Compute Lls via economic SVD decompostion
[U, S, V] = svd(Kn,0);          % Perform SVD
ss = diag(S);                   % Singular values
idxZeros = find(ss < 1E-14);    % Replace any small singular value with inf 
ss(idxZeros) = inf;
Sv = diag(1./ss);  %Inverse singular value matrix

%Least squares solution
theta = N*V*Sv*U'*Z;

%Projection matrix and residuals
P = (eye(length(Z))-U*U'); % Projection matrix
R = P*Z;                     % Residuals

%Estimate of noise and parameter variance
Cv = (R'*R)/real(trace(P));       % Noise variance
Ct = diag(N'*V*Sv*Sv*V'*N*Cv);    % Parameter variance

