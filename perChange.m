function [stats] = perChange(X,x)
% Calculates the percentage change in capacity or resisatance of battery 
%
% Mandotory input arguments
%   X:  Capacity or resitance data. Type double of size n x S. n is the
%       number of snapshot measurements. S is the number of cells
%       over which the average and standard error is calculated
%   x: Days, charge throughput or number of cycles of the snapshots. Type
%   double of size n x 1. This is used to calculate rate of degredation.
%
% Output arguments
%  stats: structure variable with fields
%         mean: mean percentage change of the samples n x 1
%         stdErr: Standard error of the percentage change, n x 1
%         confInt: 95% confidence interval, n x 1 
%         degRate: Rate of degredation for every 100 units of x
%         Xper: Percentage change matrix n x S
%
% W.D. Widanage 13/11/2018 (Grinding again)

[snapshots,nC] = size(X); % Get the size of x
S = sum(X>1,2); % Total number of valid cells
Xnorm = X.*repmat(1/mean(X(1,:),'omitnan'),snapshots,nC)*100;
x = x(:); % Vectorise x

stats.mean = mean(Xnorm,2,'omitnan');
stats.stdErr = std(Xnorm,0,2,'omitnan')./sqrt(S);
stats.confInt = stats.stdErr * 1.96; 
stats.fracUnc = stats.stdErr./stats.mean;
stats.median = median(Xnorm,2,'omitnan'); 
stats.min = min(Xnorm,[],2);
stats.max = max(Xnorm,[],2);
stats.l = stats.median - stats.min;
stats.u = stats.max - stats.median;


% Calculate rate of degredation by assuming a linear fade/increase
K = [x(:), ones(size(x))];
degRatePara = K\stats.mean;
stats.degRate = degRatePara(1)*100;    % Get the slope of degredation and scale for every 100 units
end

