function [dataSeg]  = extractStepSegments(data,mode,stepNum)

idxStepNum = mode == stepNum;


fldNames = fieldnames(data);

if idxStepNum(1) ~= 0 && idxStepNum(end) == 0
    idxStepNum = [0; idxStepNum];
    numStepNumOccStart = find(diff(idxStepNum) == 1);
    numStepNumOccEnd = find(diff(idxStepNum) == -1);
elseif idxStepNum(end) ~= 0 && idxStepNum(1) == 0
    idxStepNum = [idxStepNum; 0];
    numStepNumOccStart = find(diff(idxStepNum) == 1);
    numStepNumOccEnd = find(diff(idxStepNum) == -1);
    numStepNumOccEnd(end) = length(data.(fldNames{1}));
elseif idxStepNum(1) ~= 0 && idxStepNum(end) ~= 0
    idxStepNum = [0; idxStepNum; 0];
    numStepNumOccStart = find(diff(idxStepNum) == 1);
    numStepNumOccEnd = find(diff(idxStepNum) == -1);
    numStepNumOccEnd(end) = length(data.(fldNames{1}));
else
    numStepNumOccStart = find(diff(idxStepNum) == 1);
    numStepNumOccEnd = find(diff(idxStepNum) == -1);
end



for jj = 1:length(numStepNumOccStart)
    idxStep = [numStepNumOccStart(jj):numStepNumOccEnd(jj)];
    
    for ii = 1:length(fldNames)
        tmpVar = data.(fldNames{ii});
        dataSeg(jj,1).(fldNames{ii}) = tmpVar(idxStep);
    end
end
end