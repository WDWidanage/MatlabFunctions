function [thetaOpt] = OCVPredictSandBox(x,y,varargin)
% Predict the OCV value by fitting ym = U + sum_i a_iexp(b_ix) function,
% subject to U + sum_i a_i = y(0).
%
% Input arguments:
%   x: time vector size N x 1
%   y: Voltage data set with relaxation size N x 1
%
% Output:
%   thetaOpt: Optimum parameters based on constrained optimisation
%             theta = [U, a1;..;an, b1;...;bn] size nExp x 2
%             The estimated OCV is theta(1)
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 25/12/2016 (The Thing That Should Not be)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only

p = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
addRequired(p,'x', @isnumeric);
addRequired(p,'y', @isnumeric);

% Optional parameters
addParameter(p,'nExp',[2:3])

% Re-parse parObj
parse(p, x, y, varargin{:})

if y(end) > y(1)            % A discharge pulse
    yGP = -y + max(y);
else
    yGP = y - min(y);
end
[thetaGP, optModelOrder] = GraphPeelOpt(x,yGP,p.Results);
optNExp = optModelOrder.optModelOrder;
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',2000,'SpecifyObjectiveGradient',false);

problem.solver = 'fmincon';
objFcn = @(theta)sumSqErr(theta,x,y);
problem.objective = objFcn;
if y(end) > y(1)            % A discharge pulse
    problem.x0 = [max(y); -thetaGP(1:optNExp); thetaGP(optNExp+1:end)];
%     problem.lb  = [max(y)+1E-3;-inf(optNExp,1)];
else
    problem.x0 = [min(y); thetaGP];
end

nonLinCon =  @(theta)initalPointConst(theta,y);
problem.nonlcon = nonLinCon;

[thetaOpt, objFcnFmin]= fmincon(problem);

[thetaGA, objFcnGA] = ga(objFcn,length(problem.x0),[],[],[],[],[],[],nonLinCon);
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


function [thetaOpt0,outOpt] = GraphPeelOpt(x,y,varargin)
% Perfom graph peeling and optimisation to determine number of sum of
% exponentials in data. y = sum aexp(bx). Number of data points should be
% more than 4.
%
% Input arguments:
%   x: exponent
%   y: function value
%
% Optional input arguments. Create a structure variable with the following field:
%   nExp: Set the maximum number of exponenetials to fit up to. Default
%         nExp = [1:3]
%   plotGP: Set plotGP to 'y' to plot graph-peeling curve. Default 'n'
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
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 16/10/2016 (Calm)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only


p = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
addRequired(p,'x', @isnumeric);
addRequired(p,'y', @isnumeric);

% Optional parameters
% addParameter(p,'plotGP','n')
addParameter(p,'nExp',[1:3]);

% Re-parse parObj
parse(p, x, y, varargin{:})

x = p.Results.x(:);
y = p.Results.y(:);

nData = length(x);

nExp = p.Results.nExp;                  % Maximum number of exponentials
loopCnt = 0;                            % Initialise
for nn = nExp                           % Loop over maximum possible number of exp segments
    loopCnt = loopCnt + 1;              % Loop counter for each nExp
    z = y;                              % Make a copy of the data for graph peeling
    idxZero = find(z<=0);
    z(idxZero) = 1E-3;
    nPts = floor(nData/nn);
    
    for pp = nn:-1:1                    % Loop over each exp segment
        lnz = log(z);                   % Log of raw data
        sIdx = nPts*(pp-1) + 1;
        if pp == nn
            eIdx = nData;
        else
            eIdx = sIdx + nPts - 1;
        end
        
        lnzSeg = lnz(sIdx:eIdx);
        xSeg = x(sIdx:eIdx);
        K = [xSeg,ones(length(lnzSeg),1)];
        pLin = Lls(K,lnzSeg);             % Slope and intercept
        
        
        b(pp,1) = pLin(1);                % Exponent (slope of line)
        a(pp,1) = exp(pLin(2));           % Coefficient
        
        % Peel off exponential term
        z = z - a(pp)*exp(b(pp)*x);
        
        idxZero = find(z<=0);
        z(idxZero) = 1E-3;
    end
    
    % Evalaute CF
    ym = SumOfExp([a;b],x);
    cF(loopCnt,1) = norm(y-ym)^2;
    
    
    % Peform parameter optimastion for given graph-peel values
%     fh = @SumOfExp;
%     optionsExp.Jacobian = 'on';
%     optionsExp.termMsg = 'n';
%     optionsExp.iterMax = 100;
%     [thetaOptTmp,infoTheta] = LMAlgorithm_varIdx(fh, [a;b], x, y, optionsExp);
%     aOpt = thetaOptTmp(1:nn);
%     bOpt = thetaOptTmp(nn+1:end);
%     
%     %     options = optimset('TolX',0.1);
%     %     fh = @(lambda)sumExpErr(lambda,x,y);
%     %     bOpt = fminsearch(fh,b,options);
%     %     [sse,aOpt] = sumExpErr(bOpt,x,y);
%     
%     % Store optimised and intial values
%     theta(loopCnt,1) = {[aOpt, bOpt]};
    theta0(loopCnt,1) = {[a,b]};
%     cF(loopCnt,1) = infoTheta.cF_iter(end);
%     thetaFracErr(loopCnt,1) = {infoTheta.fracErr};
    
%     % Corrected AIC
%     dof = 2*nn+1;
%     AIC(loopCnt,1) = nData*log(infoTheta.cF_iter(end)/nData) + 2*dof + 2*dof*(dof+1)/(nData-dof-1);
%     
%     % Note any poitive exponents to avoid and flag warning
%     if any(bOpt>0)
%         outOpt.msg = sprintf('Positive exponent for nExp = %d', nn);
%         fprintf([outOpt.msg,'\n'])
%         avoidExp(loopCnt,1) = true;
%     else
%         avoidExp(loopCnt,1) = false;
%     end
    
end

% % Only select valid optimised values
% thetaValid = theta(~avoidExp,1);
% thetaFracErrValid = thetaFracErr(~avoidExp,1);
% theta0Valid = theta0(~avoidExp,1);
% nExpValid = nExp(~avoidExp);
% cFValid = cF(~avoidExp);
% AICValid = AIC(~avoidExp);


% Find optimum number of exponentials based on cost function or AIC
[~,idx] = min(cF);
optModelOrder = nExp(idx);
% thetaOpt = thetaValid{idx,1};
% thetaOptFracErr = thetaFracErrValid{idx,1};
thetaOpt0 =  theta0{idx,1};
thetaOpt0 = thetaOpt0(:);

% Largest time constant and standard deviation in optimum solution
% [minExponent,idxMin] = min(-thetaOpt(:,2));
% timeConstMax = 1./minExponent;
% fracErrMinExponent = thetaOptFracErr(optModelOrder + idxMin);
% stdMaxTimeConst = timeConstMax*fracErrMinExponent;
% stdMaxTimeCont = (-1./(minExponent+stdMinExponent)) - timeConstMax;


outOpt.theta0 = theta0;                              % Save all optimised values for each of the nExp
outOpt.cF = cF;                                      % Save the cost-function calue for each of the nExp
outOpt.optModelOrder = optModelOrder;                % Save optimum model order
% outOpt.slowestTimeConst = timeConstMax;              % Save the slowest time constant
% outOpt.slowestTimeConstStd = stdMaxTimeConst;        % Save the slowest time constant std


% 
% % Plot cost-fucntion value for each nExp and optimum fit for smallest cF
% if strcmpi(p.Results.plotGP,'y')
%     figure
%     subplot(2,1,1)
%     plot(nExpValid,AICValid,'-o')
%     xlabel('Number of exponentials'); ylabel('AIC value')
%     
%     subplot(2,1,2)
%     yOpt = SumOfExp(thetaOpt(:),x);
%     yhat = SumOfExp(thetaOpt0(:),x);
%     
%     plot(x,y,'-o',x,yhat,'- x',x,yOpt,'. -');
%     xlabel('x');
%     ylabel('y');
%     legend('Measured','Graph peel','Optimised');
%     title(sprintf('Fit for nExp = %g',nExpValid(idx)))
% end

end



function [y,J] = SumOfExp(theta,x)
% Sum of exponentials y = sum_i^M aexp(bx). i = 1...M.
%
% Input arguments:
%   theta: [a1,...,aM,b1,...,bM], size M+1 x 1 with M the number of
%          exponentials
%       x: Exponent values, size N x 1
%
% Output
%   y: function value size N x 1
%   J: Jacobian matrix size N x M
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 22/10/2016 (Black!)
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

nExp = length(theta)/2;
a = theta(1:nExp);
b = theta(nExp+1:end);


y = zeros(nData,1);

for ii = 1:nExp
    yTmp = a(ii)*exp(b(ii)*x);
    y = y + yTmp;
    J1(:,ii) = exp(b(ii)*x);
    J2(:,ii) = a(ii)*x.*exp(b(ii)*x);
end

J = [J1,J2];

end
