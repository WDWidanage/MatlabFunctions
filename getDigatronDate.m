function endDate = getDigatronDate(currentFile)

fid = fopen(currentFile);
frewind(fid);

rowStr = '';
while ~strcmp(rowStr,'End Time')
    tLine = fgetl(fid);
    rowTmp = split(tLine,',');
    rowStr = rowTmp{1};
end

endDate = datetime(rowTmp{2});
endDate.Format = 'dd-MM-yyyy';

fclose(fid);
end