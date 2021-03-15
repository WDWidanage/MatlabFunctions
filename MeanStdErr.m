function [stats] = MeanStdErr(X)
% Calculates the mean standard error and confidence interval of X. 
%
% Mandotory input arguments
%   X: Data to average. Type double of size n x S. n is the number of independant measurements. S is the number of samples
%   over which the average and standard error is calculated
%
% Output arguments
%  stats: structure variable with fields
%         mean: mean of the samples n x 1
%         stdErr: Standard error, n x 1
%         confInt: 95% confidence interval, n x 1 
%         fracUnc: Fractional uncertainity, n x 1
%         median: Sample meadian n x 1
%         min: Sample minimum
%         max: Sample maximum
%         l:  Median - minimum
%         u: Maximum - median
%
% W.D. Widanage 09/11/2018 (Grinding)

S = sum(X>1,2); % Total number of valid cells
stats.mean = mean(X,2,'omitnan');
stats.stdErr = std(X,0,2,'omitnan')./sqrt(S);
stats.confInt = stats.stdErr * 1.96; 
stats.fracUnc = stats.stdErr./stats.mean;
stats.median = median(X,2,'omitnan'); 
stats.min = min(X,[],2);
stats.max = max(X,[],2);
stats.l = stats.median - stats.min;
stats.u = stats.max - stats.median;
end

