function [theta] = Tls(K,Z,varargin)
%
% Computes a numerically stable total linear least squares estimate
% and returns the optimal estimate and its variance. Parameters will be 
% non-unique if K is rank deficient.
%
% Inputs (mandatory):
%   K: Regresor matrix, size n x m
%   Z: Output vector, size n x 1
%
% Outputs:
%   theta: Optimum total least squares parameter estimate, size m x 1
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 09-05-2020 (6LACK - Prblms)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only

Z = Z(:); % Vectorise
parObj = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
addRequired(parObj,'K', @isnumeric);
addRequired(parObj,'Z', @isnumeric);

% Re-parse parObj
parse(parObj,K,Z,varargin{:})

K = parObj.Results.K;
Z = parObj.Results.Z;

% Create augmented data and regressor matrix
Ktilda = [K,Z];                 % Augmented matrix is rank(K)+1

%Compute SVD of augmented matrix
[~, ~, V] = svd(Ktilda);  % Perform SVD

Ve = V(:,end); % Last column forms a basis vector of the null space of a rank(K)augmented matrix

% Tls solution is of the folrm [theta -1]'
Ve_norm = -Ve./(Ve(end));
theta = Ve_norm(1:end-1);














