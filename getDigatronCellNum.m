function cellNum = getDigatronCellNum(currentFile,varargin)

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments

addParameter(parObj,'numDigits',3)
parse(parObj,varargin{:});

fid = fopen(currentFile);
frewind(fid);

numDigits = parObj.Results.numDigits;

rowStr = '';
while ~strcmp(rowStr,'Circuit')
    tLine = fgetl(fid);
    rowTmp = split(tLine,',');
    rowStr = rowTmp{1};
end

cicuitStr = rowTmp{2};
cellNum = str2double(cicuitStr(end-numDigits+1:end));
fclose(fid);

end