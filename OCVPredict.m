function [thetaOpt] = OCVPredict(x,y,theta0)
% Predict the OCV value by fitting ym = U + sum_i a_iexp(b_ix) function,
% subject to U + sum_i a_i = y(0).
%
% Input arguments:
%   x: exponent
%   y: function value
%
% Optional input arguments. Create a structure variable with the following field:
%   plotFit: Set plotFit to 'y' to plot graph-peeling curve. Default 'n'
%
% Outputs:
%   thetaOpt: Overall optimum parameters based on the lowest cost-function.
%             theta = [a1;..;an, b1;...;bn] size nExp x 2
%   outOpt:  A structured varaible with the following fields as optional extras
%            theta:             Optimised set of parameters for each nExp
%            cF:                Value of cost-function for each nExp
%            modelOrder:        Optimum number of exponentials (model
%                               order)
%            slowestTimeConst:  Slowest time constant for the optimum model
%                               order
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 25/12/2016 (The Thing That Should Not be)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only




options = optimoptions('fmincon','MaxFunctionEvaluations',1000,'Algorithm','sqp','SpecifyObjectiveGradient',false);
problem.options = options;
problem.solver = 'fmincon';
problem.objective = @(theta)sumSqErr(theta,x,y);
problem.x0 = theta0;
problem.nonlcon = @(theta)initalPointConst(theta,y);

thetaOpt = fmincon(problem);
end

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


% Equatlity constraint function
function [y0InEq,y0Eq] = initalPointConst(theta,y)

y0InEq = []; % Set inequlatiy output argument to zero for fmincon function

nExp = length(theta(2:end))/2;  % Number of exponentials
U = theta(1);
a = theta(2:nExp+1);
y0Eq = U + sum(a) - y(1);       % y(0) = U + sum_i a_i
end

% Cost-function
function [sse,gradF] = sumSqErr(theta,x,y)

[ym,J] = SumOfExpOffSet(theta,x);
res = y - ym;
sse = norm(res)^2;
gradF = J'*res;
end
