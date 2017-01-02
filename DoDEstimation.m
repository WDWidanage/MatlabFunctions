function [deltaSoC, timePts] = DoDEstimation(u,time,dodThres,alpha,diag)
%
% Finds and returns all the depth-of-discharge values that occur above a
% specified threshold value.
%
% The code finds the turning points of the signal u, by using Total
% Variation Regularisation, and calculates the change in SoC between the
% identified turning points. If this value is higher than the specified
% threshold, it is stored and returned.
%
% Mandotory input arguments:
%   u:    Soc signal specified to be between 0 - 100%, size N x 1
%   time: Corresponding time signal of u, size N x 1
%
% Optional arguments:
%   dodThres: Threshold value to classify significant change in SoC,
%             default 40% as threshold, size 1 x 1
%   alpha:    Regularisation scalar factor, default 1E-3, size 1 x 1
%   diag:     Set diag to 1 for diagnosis,and to 0 for no diagnosis, default ,
%             size 1 x 1
%
% Output arguments:
%   deltaSoC: Vector of depth-of-discharge values above the specified
%             threshold, size varN x 1
%   timePts:  Vector of corresponding time intants at which SoC exceeds the
%             specified threshold, size varN x 1
%
% W. D. Widanage 16-09-2014 (think^3)

% Vectorise input and time
u = u(:);
time = time(:);

N = length(u); % Number of points

if nargin < 3 || isempty(dodThres)
    dodThres = 40;
end
if nargin < 4 || isempty(alpha)
    alpha = 1E-3;
end
if nargin < 5 || isempty(diag)
    diag = 0;
end
dt = (time(2) - time(1))*ones(N,1);

% Call derivative function to obtain smooth derivative of SoC
du = Derivative(u,alpha,dt,[],[],diag);

idxTP = find(abs([0;diff(sign(du))]) == 2); % Find the indices of the zero crossing points of the derivative

% Include the starting and end point of the soc signal along with the its
% turning points and calculate delta soc.
idx = [1; idxTP; N];

% Calculate the change in SoC or delta SoC
deltaSoC = diff(u(idx));

% Get the corresponding time points and eliminate first time instant since
% deltaSoC vector is one element shorter from the diff operation
timePts = time(idx); timePts(1) = [];

% Find delta SoC values above the specified threshold
idxAboveThres = find(abs(deltaSoC) > dodThres);
deltaSoC = deltaSoC(idxAboveThres);
timePts = timePts(idxAboveThres);

% Get indices of the timePts relative to the original time vector
[~,idxTimePts] = ismember(timePts,time);

if isempty(idxTimePts)
    warning(['No turning points found that exceed the specified threshold of ',num2str(dodThres), 'units with the regularisation parameter alpha at ',num2str(alpha)])
end

% Plot signal, simple finite difference, regularised solution and mark the
% turning points with red circles, if diagnosis is true
if diag
    figure
    subplot(3,1,1)
    plot(time, u, time(idxTP), u(idxTP), 'ro',...
        time(1), u(1), 'go', time(N),u(N),'go'); hold on;
    plot(timePts, u(idxTimePts), 'ro','MarkerFaceColor','r'); % Plot the soc points which exceed the threshold value
    ylabel('Signal to differentiate u')
    
    subplot(3,1,2)
    plot(time, [0;diff(u)]/dt);
    ylabel('Simple finite difference, dudt');
    
    subplot(3,1,3)
    plot(time, du, '. -', time(idxTP), du(idxTP), 'ro'); grid on; hold on
    plot(timePts, du(idxTimePts), 'ro','MarkerFaceColor','r'); % Plot the dsoc points which exceed the threshold value
    xlabel('Time (s)'); ylabel('TV regularised solution, dudt');
end

