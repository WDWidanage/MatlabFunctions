function containtsText = getDigatronExp(currentFile,nvp)

arguments
    currentFile {mustBeFile}
    nvp.Text {mustBeText} = ""
end

fid = fopen(currentFile);
frewind(fid);

rowStr = '';
while ~strcmp(rowStr,'Program')
    tLine = fgetl(fid);
    rowTmp = split(tLine,',');
    rowStr = rowTmp{1};
end

containtsText = contains(rowTmp{2},nvp.Text);

fclose(fid);
end