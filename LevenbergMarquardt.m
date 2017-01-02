function [theta,results] = LevenbergMarquardt(fh, theta0, u, d, varargin)
%
% Perform nonlinear least squares optimisation using the
% Levenberg-Marquardt method
%
% Minimised cost-function f(theta)
%   f(theta) = sum ((G(theta)_i-d_i)/s_i)^2  i = 1...M
%
% Mandotory input argumetns
%   fh: function handle of nonlinear model. fh = @(theta, u)fcn(theta,u,optArg1,...,optArgN).
%       Nonlinear model function should return two outputs, the model
%       output and Jacobian as second (Jacobian is optional)
%   theta0: Initial starting point of model parameters, size N x 1
%   u: Input signal to simulate nonlinear model, size M x 1
%   d: Measured output data, size M x 1
%
% Optional arguments. Create a structure variable with the following fields:
%   Jacobian: Specify as 'on' if model function returns the Jacobian or to 'off' to approximate by finte forward difference, default Jacobian ='off'.
%   s: Residual weights, normally std of noise, default s = 1, size M x 1
%   iterMax: Maximum number of iterations, default iterMax = 1000, double size 1 x 1
%   TolDT:  Termination tolerance of parameter update, default TolDT = 1E-6, double size 1 x 1
%   diagnosis: Set daignosis to 'on' to plot cost-function and lambda vs iterations
%
% Output arguments:
%   theta: Optimised parameter vector, size N x 1
%   results: Sturcture variable with
%            results.cF_iter: cost-function value at each iteration 
%            results.L_iter: Lambda weight at each iteration 
%            results.Msg: Reason for iteration termination
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 14/10/2015 (Highway to hell!!)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only


p = inputParser; % Create an input parse object to handle positional and property-value arguments
theta0 = theta0(:);
nData = length(u);
nPara = length(theta0);

% Create variable names and assign default values after checking the value
addRequired(p,'fh', @checkFunctionHandle);
addRequired(p,'theta0', @isnumeric);
addRequired(p,'u', @isnumeric);
addRequired(p,'d', @isnumeric);

% Optional parameteres
addParameter(p,'Jacobian','off')
addParameter(p,'s',ones(size(u)))
addParameter(p,'iterMax',1000,@isnumeric)
addParameter(p,'TolDT',1E-6,@isnumeric)
addParameter(p,'TolDCF',1E-8,@isnumeric)
addParameter(p,'diagnosis','off')

% Re-parse parObj
parse(p,fh, theta0, u, d, varargin{:})


% Initialise
theta_prev = p.Results.theta0;
% Evalaute model function and Jacobian for initial parameter values
if ismember(p.Results.Jacobian,{'On','on'})
    [y,J_prev] = fh(theta_prev,p.Results.u);
end
[cF_prev,F_prev] = costFunctionEval(y,p.Results.d,p.Results.s);

cF = cF_prev*10; % Induce that the present cost-function is worse than previous one with the assumed value of lambda
iterUpdate = 1;
deltaT = 1E6;
deltaCF = abs(cF - cF_prev);
lambda = 10;
innerLoop = 1; % Used for de-bugging purposes

cF_iter(iterUpdate,1) = cF;
L_iter(iterUpdate,1) = lambda;

% Start solution update
while norm(deltaT) > p.Results.TolDT && iterUpdate <= p.Results.iterMax  %&& deltaCF > p.Results.TolDCF
    while cF > cF_prev % Increase lambda and re-evaluate parameter update
        lambda = lambda*10;
        deltaT = parameterUpdate(J_prev,F_prev,lambda);     % Calucuate parameter update
        theta = theta_prev + deltaT;                        % Update parameter estimate
        
        % Evalaute model function and Jacobian for updated parameter
        if ismember(p.Results.Jacobian,{'On','on'})
            [y, J] = fh(theta,p.Results.u);
        end
        [cF, F] = costFunctionEval(y,p.Results.d,p.Results.s);
        innerLoop = innerLoop + 1;
    end
    
    cF_prev = cF;
    theta_prev = theta;
    J_prev = J;
    F_prev = F;
    
    lambda = lambda/10;
    deltaT = parameterUpdate(J,F,lambda);       % Calucuate parameter update
    theta = theta + deltaT;                     % Update parameter estimate
    
    % Evalaute model function and Jacobian for updated parameter
    if ismember(p.Results.Jacobian,{'On','on'})
        [y,J] = fh(theta,p.Results.u);
    else
        [y] = fh(theta,p.Results.u);
    end
    [cF, F] = costFunctionEval(y,p.Results.d,p.Results.s);
    deltaCF = abs(cF - cF_prev);
    
    cF_iter(iterUpdate,1) = cF;
    L_iter(iterUpdate,1) = lambda;    
    iterUpdate = iterUpdate + 1;
    
end % End of main while iterative loop

iterUpdate = iterUpdate - 1; % Reduce iteration count by one when loop is exited

sCF = cF^2/(nData-nPara);
% covTheta = (J'*J)\eye(nPara,nPara)*sCF;
covTheta = CovTheta(sCF,J);    % Parameter variance

if ismember(p.Results.diagnosis,{'On','on'})
    diagnosis(cF_iter,L_iter,iterUpdate)
end
% idx = [norm(deltaT) < p.Results.TolDT, iterUpdate == p.Results.iterMax, deltaCF < p.Results.TolDCF];
% termStr = {[' Parameter update is smaller than specified tolerance, TolDT = ', num2str(p.Results.TolDT),'.'],...
%            [' Maximum iteration reached, iterMax = ', num2str(p.Results.iterMax),'.'],...
%            [' Change in cost-function lower than TolDCF = ',num2str(p.Results.TolDCF),'.']};

idx = [norm(deltaT) < p.Results.TolDT, iterUpdate == p.Results.iterMax];
termStr = {[' Parameter update is smaller than specified tolerance, TolDT = ', num2str(p.Results.TolDT),'.'],...
           [' Maximum iteration reached, iterMax = ', num2str(p.Results.iterMax),'.']};
disp(['Iteration terminated: ',termStr{idx}]);

results.covTheta = covTheta;
results.stdTheta = sqrt(diag(covTheta));
results.cF_iter = cF_iter;
results.L_iter = L_iter;
results.Msg = ['Iteration terminated: ',termStr{idx}];
end

function valid = checkFunctionHandle(fh)
testFH = functions(fh);
if testFH.function
    valid = true;
else
    valid = false;
end
end


function [cF,F] = costFunctionEval(y,d,s)
F = (y-d)./s;       % Weighted residual
cF = norm(F)^2;     % Cost-function
end


function deltaT = parameterUpdate(J,F,lambda)
K = ((J')*J + lambda*diag(diag((J')*J)));       % Create LM regressor matrix (cost-function Hessian + Steepest descent)
Z = (-J')*F;                                    % Negative cost-function gradient
deltaT = Lls(K,Z);                          %  all numerically stable linear least squares method
end


function diagnosis(cF,L,I)
figure()
semilogx([0:I-1],cF,'.-')
xlabel('Iteration number')
ylabel('Cost-fucntion (y-d/s)^2')

figure()
plot([0:I-1],L,'.-')
xlabel('Iteration number')
ylabel('Steepest descent lambda factor')
end



function theta = Lls(K,Z)
% Computes a numerically stable linear least squares estimate
% and returns the optimal estimate, its variance and an estimate of the
% noise variance
%
% Inputs (mandatory): 
%   K: Regresor matrix, size n x m
%   Z: Output vector, size n x 1
%
% Outputs:
%   theta: Optimum parameter estimate, size m x 1


% Normalise K with the 1/norm of each column
Knorm = sqrt(sum(abs(K).^2));
idxZeros = Knorm<1E-14;
Knorm(idxZeros) = 1;
N = diag(1./Knorm);
Kn = K*N;

% Compute Lls via economic SVD decompostion
[U, S, V] = svd(Kn,0);          % Perform SVD
ss = diag(S);                   % Singular values
idxZeros = ss < 1E-14;    % Replace any small singular value with inf 
ss(idxZeros) = inf;
Sv = diag(1./ss);  %Inverse singular value matrix

% Least squares solution
theta = N*V*Sv*(U')*Z;

end

function ct = CovTheta(s,J)
% Numerically stable calculation of the parameter covaraince matrix.
%
% Inputs
%   s: Sum of squared residuals/(nDataPts - nPara
%   J: Jacobian matrix
%
% Outputs
%   ct: Parameter covariance matrix

% Normalise J with the 1/norm of each column
Jnorm = sqrt(sum(abs(J).^2));
idxZeros = Jnorm<1E-14;
Jnorm(idxZeros) = 1;
N = diag(1./Jnorm);
Jn = J*N;

[~, S, V] = svd(Jn,0);          % Perform SVD
ss = diag(S);                   % Singular values
idxZeros = ss < 1E-14;    % Replace any small singular value with inf 
ss(idxZeros) = inf;
Sv = diag(1./ss);  %Inverse singular value matrix

ct = (N')*V*(Sv)*Sv*(V')*N*s;    % Parameter covariance matrix
end