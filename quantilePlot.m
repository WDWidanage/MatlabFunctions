function varargout = quantilePlot(X,varargin)
% Plot a quantile-quantile plot of the data set X to test for normality of
% the data.
%
% Mandotory inputs
%   X: Cell array. Each entry of X is the data set to test for normality.
%      cell array of size 1 X M. Each element is an array of varying size,
%      X_i x 1
%
%
% Optional input: Create a structure variable with the following field:
%   disp: Set disp to 'off' to avoid plotting. If 'off' then read the
%   optional ouput and look at the corelation coefficient value > 0.8 for
%   nomality.
%
% Optional output
%  varargout: A structure variable with the results of qq analaysis with
%            fields
%            xSort:     Sorted (ascending) measured data set
%            idxSort:   Indices of sorted measured datapoints
%            z:         Theoretical data points
%            corrCoeff: Sample correlation coefficient between measured and
%                       theoretical data
%
% W.D . Widanage 04-03-2016 (Deal with the Four Horsemen!!)

parObj = inputParser;

% Create variable names and assign default values after checking the value
addRequired(parObj,'x');

% Optional parameteres
addParameter(parObj,'disp','on');

% Re-parse parObj
parse(parObj,X,varargin{:})


for xx = 1:length(X)
    x = X{xx};
    nX = length(x);     % Number of samples
    mX = mean(x);       % Sample mean
    stdX = std(x);      % Sample std
    
    q = ([1:nX]-0.5)/nX;    % Calculate quantiles
    zTmp = norminv(q,0,1);  % Calculate theoretical values based on inverse standard normal distribution
    z = zTmp'*stdX + mX;     % Translate theoretical data with sample mean and std
    [xSort,idxSort] = sort(x);        % Sort data for plotting and correlation of coefficient calculation
    
    corrCoeff = (mean(xSort.*z)-mean(xSort)*mean(z))/(std(xSort)*std(z)); % Correlation of coefficient
    
    titleMsg = sprintf('Correlation coefficient r: %.2f%',corrCoeff);
    if strcmpi(parObj.Results.disp,'on')
        figure()
        plot(z,z,'r',z,xSort,'x')
        xlabel('Expected data'); ylabel('Measured data')
        title(titleMsg)
    end
    
    qqResults.(['DataSet',num2str(xx)]).xSort = xSort;
    qqResults.(['DataSet',num2str(xx)]).z = z;
    qqResults.(['DataSet',num2str(xx)]).corrCoeff = corrCoeff;
    qqResults.(['DataSet',num2str(xx)]).idxxSort = idxSort;
    
    
end

varargout{1} = qqResults;

end

