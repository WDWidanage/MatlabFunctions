function [idx_Exc,idx_ExcRelax] = find_exc_segments(current,varargin)
% This function is for cell pulse testing. Based on a current load profile the
% function returns the index pairs of the excitation segments and
% excitation with relaxation segements. The function provides all the
% excitation and relaxation segments witin the first and last rest
% intervals
%
% Mandotory input arguments:
%   current: Current vector (A), vector size N x 1
%
% Optional input arguments. Create a structure variable with the following fields:
%   tol: tolorance level to define zero level, default 0.1A, size 1 x 1
%   plotSeg: set to 1 to plot ansd display signal segments, default 0, size 1 x 1
%   timeVec: time column of measured current. This is used for plotting, vector size N x 1
%
% Output arguments:
%   idx_Exc: Index pairs for excitation segment, matrix p x 2
%   idx_ExcRelax: Index pairs for excitation segment with relaxation, matrix p x 2
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 24/09/2013 (Welcome to the jungle!!)
% All Rights Reserved

pObj = inputParser; % Create an input parse object to handle positional and property-value arguments


% Create variable names and assign default values after checking the value
addRequired(pObj,'current', @isnumeric);

% Optional parameters
addParameter(pObj,'tol',0.1,@isnumeric);
addParameter(pObj,'plotSeg',0,@isnumeric);
addParameter(pObj,'timeVec',[],@isnumeric);


% Re-parse parObj
parse(pObj,current,varargin{:})

current = pObj.Results.current;
tol = pObj.Results.tol;
plotSeg = pObj.Results.plotSeg;
timeVec = pObj.Results.timeVec;

logical_segments = -tol < current & current < tol; % logic 1 => zero segments, logic 0 => excitation segments

% Initialise
if logical_segments(1) == true
    s(1) = 1;       % initialise start index to 1 and set start vector index to 1
    jj = 1;
else
    jj = 0;        % else set start vector index to 0 since the value will arrive later in loop
end

kk = 0;            % set end vector index to 0 since the value will arrive later in loop
for ii = 2:length(current)
    logical_prev = logical_segments(ii-1);
    if logical_segments(ii) ~= logical_prev    % Detect a change in event
        if logical_segments(ii) == true        % True indicates start of zero segment
            jj = jj+1;
            s(jj) = ii;                        % Save index of start segment
        else                                   % False indicates end of zero segment
            kk = kk+1;
            e(kk) = ii-1;                      % Save index of end segment
        end
    end
end

% If last segment is a zero segmment there wont be a change in event in the for loop
% and manually set the last element of the end vector to the signal length  
if logical_segments(end) == true
    e(kk+1) = length(current);       
end

idx_Exc = [e(1:end-1)',s(2:end)'-1];
idx_ExcRelax = [e(1:end-1)',e(2:end)'];

[lEx, ~] = size(idx_Exc);
[lExR, ~] = size(idx_ExcRelax);

if  isempty(timeVec)
    timeVec = [1:length(current)]';
end

 if plotSeg
     figure
     plot(timeVec',current,'. -', timeVec(idx_Exc(:,1)),zeros(lEx,1),'g o','MarkerFaceColor','g');
     hold on; 
     plot(timeVec(idx_Exc(:,2)),current(idx_Exc(:,2)),'r o','MarkerFaceColor','r');
     plot(timeVec(idx_ExcRelax(:,2)),zeros(lExR,1),'k o');
     xlabel('Time'); ylabel('Curent')
     hold off;
 end
 

end