function [dataSeg]  = extractPulseStepSegments(data,mode,stepNum)

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


cntr = 0;
for jj = 1:length(numStepNumOccStart)
    if mode(numStepNumOccEnd(jj)+1) == stepNum + 1
        cntr = cntr + 1;
        modeTrunc = mode(numStepNumOccEnd(jj)+1:end);
        idxRestNum = modeTrunc == stepNum + 1;
        numStepNumOccEndR = find(diff(idxRestNum) == -1,1) + numStepNumOccEnd(jj);

        idxStep = [numStepNumOccStart(jj):numStepNumOccEndR];

        for ii = 1:length(fldNames)
            tmpVar = data.(fldNames{ii});
            dataSeg(cntr,1).(fldNames{ii}) = tmpVar(idxStep);
        end
    else
        continue
    end
end
end