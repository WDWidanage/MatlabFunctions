function  SegmentSignal(startEndTimes,cellFileName,dataMatrix)
% Split signals in to segments based on the start and end times given and
% save the segments in a mat file with the specifed file name
%
% Input arguments:
%   startEndTimes: Matrix of start and end absolute time in seconds for
%                  n segments Size n x 2 matrix.
%   cellFileName:  Cell variable with file names. Size  n x 1 cell
%   dataMatrix:    A matrix of concatonated time in secoinds and other
%                  signals to be segmented. Size m x s
%
% W. D. Widanage 20/02/2014 (mundane)



% Get the number of segments
[nSeg, ~] = size(startEndTimes);

% Get the time vector from the dataMatrix
timeVecTmp = dataMatrix(:,1);          % This should be in seconds

for ii = 1:nSeg
    % Get the start and end absolute times in seconds
    seTime = startEndTimes(ii,:);
    idxTmp = find(timeVecTmp <= seTime(1));
    idxStart = idxTmp(end);
    idxTmp = find(timeVecTmp >= seTime(2));
    idxEnd = idxTmp(1);
    
    matFileName = cellFileName{ii};
    timeStart = dataMatrix(idxStart,1);                     % Start time
    timeVec = dataMatrix([idxStart:idxEnd],1) - timeStart;  % Relative Time vector (s)
    currentVec = dataMatrix([idxStart:idxEnd],2);           % Current (A)
    volVec = dataMatrix([idxStart:idxEnd],3);               % Voltage (V)
    powVec = dataMatrix([idxStart:idxEnd],4);               % Power (W)
    capVec = dataMatrix([idxStart:idxEnd],5);               % Capacity (Ah)
    eneVec = dataMatrix([idxStart:idxEnd],6);               % Energy (Wh)
    tempVec = dataMatrix([idxStart:idxEnd],7);              % Temperature (degC)
    
    % Save data in a mat file
    save(matFileName,'timeVec','currentVec','volVec','powVec','capVec','eneVec','tempVec')
    
end
end

