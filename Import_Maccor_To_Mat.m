function Import_Maccor_To_Mat(varargin)
% Syntax of use
%   Import_Maccor_To_Mat()
%   Import_Maccor_To_Mat(pwd)
%   Import_Maccor_To_Mat(pwd,pwd)
%   Import_Maccor_To_Mat(pwd,pwd,[headerRowPosition,dataRowPosition])
%   Import_Maccor_To_Mat(pwd,pwd,[headerRowPosition,dataRowPosition],'saveNames',{matFileName1,...,matFileNameN})
%   Import_Maccor_To_Mat('currSgn',1)
%   Import_Maccor_To_Mat('saveNames',{matFileName1,...,matFileNameN})
%
%
% Description of function
%   Import data from Bitrode CSV file(s) and save as MAT file(s).
%   If there is no input, or the first input is an empty array, a file
%   browser will open so that files can be selected.
%
%   Input 1 (optional) Can be any of the following:
%       An empty array (see above)
%       A single file name as a string;
%       A cell array of file names;
%       A single folder name as a string.
%       Default is []
%
%   Input 2 (optional): a string for a folder to save the file(s) in. If an
%   empty array is passed in, a GUI will open to select a folder. Default
%   is [].
%
%   Input 3 (optional): a 2 element numeric vector containing the row number
%   which has the column names, and the row at which the main data begins.
%   Common examples include [14,16], [1,3] and [0,1]. If the first element is
%   set to 0 then no column names are assumed. Default is [].
%
%   Input 4 (optional): A property-value input value of 1 or -1.
%       Property: currSgn
%       Value: Depending on how Maccor export settings are set, this property should be set to
%   1, if dicharge current is assumed positive when exporting the CSV file or to -1 if
%   exported assuming negative for discharge. Default is 1.
%
%   Input 5 (optional): A property-value input specifying file names to use
%   when saving the converted csv file.
%       Property: saveNames
%       Values: A single file name as a string or cell array of file names
%   If no filenames are specified in Input 1, saveNames will be used to save results alphabetically
%   If filenames are specified, saveNames will correspond to order of
%   filenames in Input 1. Default is [].
%
% The file name to convert can contain no path or a partial path. File names don't
% have to have the .CSV extension. If no save path is passed to the
% function, a file browser will open to let the user select a folder.
% Currently the import function can only handle the 15 row header
% spreadsheet option from the Bitrode. Don't edit the CSV files in Excel
% because it saves it as a different type of (compressed) CSV and the
% header changes so the code errors out.
%
%
% Original script by: Thomas Bruen 07-May-2015 11:16:25
% Edited to handle input arguments and saving names by: W.D. Widanage 09-July-2015
% Edited to handle csv files with problamatic headers and time columns by W. D. Widanage 21-Dec-2015
% Edited to handel Maccor txt file W. D. Widanage 23-Mar-2019
% Edited to handel Maccor meta info  W. D. Widanage 26-Oct-2019

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
% validity - 09-July-2015
addOptional(parObj,'csvTxtDetails',[], @checkCSVTxtDetails);
addOptional(parObj,'savePath',[], @checkSavePath);
addOptional(parObj,'hdPos',[], @checkHdPos);
addParameter(parObj,'currSgn',1, @isnumeric);
addParameter(parObj,'saveNames',[], @checkSaveNames);
addParameter(parObj,'fileType','.csv'); % Assume default Maccor file type is csv
addParameter(parObj,'metaPos',[], @checkMetaPos);


% Re-parse parObj - 09-July-2015
parse(parObj,varargin{:});


if isempty(parObj.Results.csvTxtDetails)
    %open GUI to find files
    [fileNames,filePath,~] = uigetfile('*.csv','Select CSV File(s)','MultiSelect','On');
    if isnumeric(fileNames) %when cancelled, UIGETFILE returns 0
        warning('No files selected')
        return
    elseif ischar(fileNames)
        fullFileNames = {[filePath,fileNames]};
    else
        fullFileNames = cellfun(@(f) sprintf('%s%s',filePath,f),fileNames,'UniformOutput',false)';
    end
else
    % Could be a single file, or a cell array of files or directory. Pass to
    % function to process and return a cell array of file names with path
    fullFileNames = Input_Process(parObj);
end

%check all files are CSVs and exist.
fullFileNames = Check_Filenames(fullFileNames,parObj);

%specify a directory to save MAT files in:
if  ~isempty(parObj.Results.savePath)
    savePath = parObj.Results.savePath;
    if ~ischar(savePath) || ~isdir(savePath)
        error('Save path cannot be found')
    end
else
    %no location requested: open file browser to select a folder
    savePath = uigetdir;
    if isnumeric(savePath) %user cancelled
        warning('No save path selected')
        return
    end
end


%loop through files and send to function to get data, then save
nFiles = numel(fullFileNames);
for n = 1:nFiles
    currentFile = fullFileNames{n};
    
    metaInfo = metaDataStr(currentFile,parObj);
    hdrData =  HdrDataRow(currentFile,parObj); % Call function determine header and data row
    
    tic
    data = Import_Data(currentFile,hdrData(1),hdrData(2),parObj);
    data.metaInfo = metaInfo;
    toc
    
    if isempty(parObj.Results.saveNames)
        [~,f,~] = fileparts(currentFile);
    else
        if ischar(parObj.Results.saveNames)
            [numSaveFiles,~] = size(parObj.Results.saveNames);
            if numSaveFiles ~= nFiles
                error('Number of csv files is not equal to number of save file names');
            end
            f = parObj.Results.saveNames;
        end
        if iscell(parObj.Results.saveNames)
            numSaveFiles = numel(parObj.Results.saveNames);
            if numSaveFiles ~= nFiles
                error('Number of csv files is not equal to number of save file names');
            end
            f = parObj.Results.saveNames{n};
        end
    end
    saveFile = [savePath,filesep,f,'.mat'];
    save(saveFile,'-struct','data')
end
%%%%%% End of main function %%%%%%%




% Check validty of inputs - 09-July-2015
function valid = checkCSVTxtDetails(csvTxtDetails)
% Check if csvTxtDetails is a string or cell array of strings
if strcmp(csvTxtDetails,'saveNames')
    valid = false;
elseif ischar(csvTxtDetails)
    valid = true;
elseif all(cellfun(@ischar,csvTxtDetails))
    valid = true;
else
    error('Argument must be a string, or cell array of strings');
    valid = false;
end




% Check validty of inputs - 09-July-2015
function valid = checkSavePath(savePath)
% Check if savePath is a valid string and directory
if strcmp(savePath,'saveNames')
    valid = false;
elseif ischar(savePath) || isdir(savePath)
    valid = true;
else
    error('Path cannot be found')
    valid = false;
end




% Check validty of inputs - 09-July-2015
function valid = checkHdPos(hdPos)
% Check if hdPos is a valid numeric vector
if strcmp(hdPos,'saveNames')
    valid = false;
elseif isnumeric(hdPos) && numel(hdPos) == 2
    valid = true;
else
    error('Invalid values for header and data position, must be a numeric vector with two elements');
    valid = false;
end




% Check validty of inputs - 09-July-2015
function valid = checkSaveNames(saveNames)
% Check if saveNames are strings
if ischar(saveNames)
    valid = true;
elseif all(cellfun(@ischar,saveNames))
    valid = true;
else
    error('Invalid save names')
    valid = false;
end

% Check validty of inputs - 09-July-2015
function valid = checkMetaPos(metaPos)
% Check if metaPos is a valid numeric vector
if strcmp(metaPos,'metaPos')
    valid = false;
elseif isnumeric(metaPos) && numel(metaPos) == 2
    valid = true;
else
    error('Invalid values for meta data position, must be a numeric vector with two elements');
    valid = false;
end


function fullFileNames = Input_Process(arg)
% parse first input to main function.
if ischar(arg.Results.csvTxtDetails)
    argCell = {arg.Results.csvTxtDetails}; %cell so it can be treated the same as if a cell array
elseif iscell(arg.Results.csvTxtDetails)
    argCell = arg.Results.csvTxtDetails;
else
    error('Argument must be a string, or cell array of strings')
end

%get file(s) name with path for three cases: a cell array of file names, a
%single file name or a directory (in which case find the files in the
%directory).
if numel(argCell) > 1 %cell array: should be filenames
    fullFileNames = argCell;
else
    %determine if input is a directory or a single file
    if isdir(argCell{1})
        %is a directory: find all CSV files
        filePath = argCell{1};
        if strcmp(arg.Results.fileType,'.csv')
            D = dir(sprintf('%s%s%s',filePath,filesep,'*.csv'));
            fileNames = {D.name}';
            if isempty(fileNames)
                error('No Maccor csv files found in directory')
            else
                fullFileNames = cellfun(@(f) sprintf('%s%s%s',filePath,filesep,f),fileNames,'UniformOutput',false);
            end
        elseif strcmp(arg.Results.fileType,'.txt')
            D = dir(sprintf('%s%s%s',filePath,filesep,'*.txt'));
            fileNames = {D.name}';
            if isempty(fileNames)
                error('No Maccor txt files found in directory')
            else
                fullFileNames = cellfun(@(f) sprintf('%s%s%s',filePath,filesep,f),fileNames,'UniformOutput',false);
            end
        end
        
    else
        %assume it is a file with path for now (check later)
        fullFileNames = argCell(1); %keep as cell array
    end
end





function fullFileNames = Check_Filenames(fullFileNames, parObj)
%check that the file names are CSVs and all exist.
for n = 1:numel(fullFileNames)
    [p,f,e] = fileparts(fullFileNames{n});
    
    %add CSV extension if there is none
    if isempty(e) & strcmp(parObj.Results.fileType,'.csv')
        e = '.csv';
    elseif isempty(e) & strcmp(parObj.Results.fileType,'.txt')
        e = '.txt';
    end
    
    %check extension is correct
    if ~strcmp(lower(e),'.csv') & strcmp(parObj.Results.fileType,'.csv')
        error('%s is not a csv file.',fullFileNames{n});
    end
    if ~strcmp(lower(e),'.txt') & strcmp(parObj.Results.fileType,'.txt')
        error('%s is not a txt file.',fullFileNames{n});
    end
    
    %if there is no path, assume current directory
    if isempty(p)
        p = pwd; %current folder
    end
    
    %rebuild filename
    fullFileNames{n} = [p,filesep,f,e];
end

%check files all exist:
fileChk = cellfun(@(f) exist(f,'file') == 2,fullFileNames);
nChk = numel(fileChk(fileChk == 0));
missingFiles = fullFileNames(~fileChk);
if nChk
    errorFormatString = repmat('%s\n',1,nChk);
    error(sprintf('The following files cannot be found: \n%s',errorFormatString),missingFiles{:});
end



function hdrDataVec = HdrDataRow(currentFile,parObj)
% Automatically determine the start of the header and data row. Three
% possible otptions.
%  1.   If first row has the word 'Test', the bitrode csv file has the extra file
%       info. Therefore header and data rows are 14 and 16
%  2.  If first row has the word 'Current', the bitrode csv file starts with the
%      header. Therefore header and data rows are 1 and 3
%  3.  Otherwise assume first row is the data row. Therefore header and data
%      rows are 0 and 3


%open file
fid = fopen(currentFile);
frewind(fid); %just in case

fstRowStr = fgetl(fid);


if  numel(parObj.Results.hdPos) == 2
    hdrDataVec= parObj.Results.hdPos;
elseif regexp(fstRowStr,'Test')
    hdrDataVec = [14,16];
elseif regexp(fstRowStr,'Current')
    hdrDataVec = [1,3];
else
    hdrDataVec = [0,3];
end






function dataStruct = Import_Data(fileName,rowColumnNames,rowDataStart,parObj)
%main routine to import data. File is loaded using FREAD and the header is
%processed to get the column names (using a subfunction). These column
%names are used as structure field names for the data. TEXTSCAN is used to
%read in the data.

%open file
fid = fopen(fileName);
frewind(fid); %just in case

%get row of column names, if there are any. If not, handle later
if rowColumnNames > 0
    %skip over headers:
    for n = 1:rowColumnNames-1 %-1: fgetl goes to the next line
        [~] = fgetl(fid);
    end
    columnNameText = fgetl(fid); % channel names
else
    columnNameText = [];
end

%skip to data
for n = 1:(rowDataStart - rowColumnNames -1)
    [~] = fgetl(fid);
end

%get first row of data to test the format, then go back to start of data
dataStartPosition = ftell(fid);
dataFirstLine = fgetl(fid);
fseek(fid,dataStartPosition,'bof');


%get preallocated structure and channel types (string or numeric)
[dataStruct,formatString] = Process_Channel_Names(columnNameText,dataFirstLine,parObj);
nChannels = numel(dataStruct.channelNames);
dataStruct.csvInfo = dir(fileName);

if strcmp(parObj.Results.fileType,'.csv')
    data = textscan(fid,formatString,'delimiter',',');
elseif strcmp(parObj.Results.fileType,'.txt')
    data = textscan(fid,formatString,'delimiter','\t');
end
fclose(fid);



% save columns of data to structure, depending on if string or numeric
dataNames = fieldnames(dataStruct.data);
for n = 1:nChannels
    dataStruct.data.(dataNames{n}) = data{n};
    if  strcmp(dataNames{n},'State')
        stateIdx = n;
    end
    if  strcmp(dataNames{n},'Md')
        stateIdx = n;
    end      
end

% Make discharge current -ve
if( parObj.Results.currSgn == -1)
    if isfield(dataStruct.data,'Amps')
        currentTmp = dataStruct.data.Amps;
        
        isDischarge = cellfun(@strcmp,data{stateIdx},repmat({'D'},length(data{stateIdx}),1));
        dataStruct.data.Amps(isDischarge) = currentTmp(isDischarge)*parObj.Results.currSgn; %
        dataStruct.nSamples = size(dataStruct.data.(dataNames{n}),1);
    end
    
    if isfield(dataStruct.data,'Current_A_')
        currentTmp = dataStruct.data.Current_A_;
        
        isDischarge = cellfun(@strcmp,data{stateIdx},repmat({'D'},length(data{stateIdx}),1));
        dataStruct.data.Current_A_(isDischarge) = currentTmp(isDischarge)*parObj.Results.currSgn; %
        dataStruct.nSamples = size(dataStruct.data.(dataNames{n}),1);
    end
end

if isfield(dataStruct.data,'TestTime')
    if iscell(dataStruct.data.TestTime(1))
        dataStruct.data.TestTime = Time_String_To_Seconds(dataStruct.data.TestTime);
    end
end


function [dataStruct,formatString] = Process_Channel_Names(channelString,dataString,parObj)
%get channel names, and test to determine which channels are text and whcih
%are numeric.
%channelString (char array): row of spreadsheet containing channel names.
%dataString (char array): row of spreadsheet containing a row of data

% %use first line of data to see how many channels there are:
% nChannels = numel(strfind(dataString,','))+1;

%read all data channels as a string, then see which can be converted to
%numbers:
if strcmp(parObj.Results.fileType,'.txt')
    idxSpace = double(dataString) == 9; % Double value of a tab
    dataString(idxSpace) = ',';
    idxSpace = double(channelString) == 9; % Double value of a tab
    channelString(idxSpace) = ',';
end
dataTextTmp = regexp(dataString,',','split');
empData = cellfun(@isempty,dataTextTmp);
dataText = dataTextTmp(~empData);
nChannels = numel(dataText);
dataConv = cellfun(@str2double,dataText);
formatString = repmat('%f',1,nChannels);
formatString(find(isnan(dataConv))*2) = 's';

%if there aren't any column names row, then create a dummy one
if isempty(channelString)           % No headers are specified. Use default, Channel_1...Channel_N
    channelNamesTmp = sprintf('Channel_%d',1:nChannels);
    channelNames = regexp(channelNamesTmp,'Channel_[\d]{0,}','match'); % Get channel hdr names
    units = repmat({'None'},nChannels,1);
else
    
    expHdr = {'Rec#','Cyc#','Step','TestTime(s)','StepTime(s)','Amp-hr','Watt-hr','Amps',...
        'Volts','State','ES','DPt Time','ACR','DCIR','Temp 1','Temp 2','VARx1','VARx2','VARx3',...
        'VARx4','VARx5','VARx6','VARx7','VARx8','VARx9','VARx10','VARx11','VARx12','VARx13',...
        'VARx14','VARx15','FLGx1','FLGx2','FLGx3','FLGx4','FLGx5','FLGx6','FLGx7','FLGx8','FLGx9',...
        'FLGx10','FLGx11','FLGx12','FLGx13','FLGx14','FLGx15','FLGx16','FLGx17','FLGx18',...
        'FLGx19','FLGx20','FLGx21','FLGx22','FLGx23','FLGx24','FLGx25','FLGx26','FLGx27','FLGx28',...
        'FLGx29','FLGx30','FLGx31','FLGx32','EV Temp','EV Hum'};
    
    channelNamesTmp = split(channelString,',');
    channelNames = channelNamesTmp(1:nChannels)';
    
    units = repmat({'None'},nChannels,1);
       
end

%create a structure to store the data in:
structNames = matlab.lang.makeValidName(channelNames); %make channel names valid field names
uniChnlNames = UniqueChannelNames(structNames);        % Handle repeating channel names
structDummy = [uniChnlNames;cell(1,nChannels)];
dataStruct.data = struct(structDummy{:});

%add channel info to structure
dataStruct.channelNames = channelNames;
dataStruct.channelUnits = units;

function chnlHdrs = UniqueChannelNames(chnlHdrs)
nEl = length(chnlHdrs);
[~,uniIdx] = unique(chnlHdrs);
repIdx = setdiff([1:nEl],uniIdx);
for ii = 1:length(repIdx)
    chnlHdrs(repIdx(ii)) = {[chnlHdrs{repIdx(ii)},sprintf('_%d',ii)]};
end




        

function timeVec = Time_String_To_Seconds(timeCell)
N = length(timeCell);
% numDays = str2num(cellfun(@extractBefore,timeCell,repmat({'d'},N,1)));
numDaysTmp = cellfun(@extractBefore,timeCell,repmat({'d'},N,1),'UniformOutput',false);
numDays =  cellfun(@str2num,numDaysTmp);
hmsCell = cellfun(@extractAfter,timeCell,repmat({'d '},N,1),'UniformOutput',false);
durCell = cellfun(@duration,hmsCell);
timeVecSec = seconds(durCell);

% if idxZ(1) == 1
%     secR = cumsum(secTmp(idxZ(2:end)-1));
%     timeVecSec(1:idxZ(2)-1,1) = secTmp(1:idxZ(2)-1);
%     for jj = 1:length(secR)-1
%         idxTmp = idxZ(jj+1):idxZ(jj+2)-1;
%         timeVecSec(idxTmp,1) = secTmp(idxTmp) + secR(jj);
%     end
% else
%     secR = cumsum(secTmp(idxZ-1));
%     timeVecSec(1:idxZ(1)-1,1) = secTmp(1:idxZ(1)-1);
%     for jj = 1:length(secR)-1
%         idxTmp = idxZ(jj):idxZ(jj+1)-1;
%         timeVecSec(idxTmp,1) = secTmp(idxTmp) + secR(jj);
%     end
% end

timeVec = timeVecSec - timeVecSec(1) + numDays*86400;



function metaInfo = metaDataStr(currentFile,parObj)
% Extract the meta data infomation

fid = fopen(currentFile);
frewind(fid);

if numel(parObj.Results.metaPos) == 2
    metaDataSrt = parObj.Results.metaPos(1);
    metaDataEnd = parObj.Results.metaPos(2);   
    
    %skip to meta data
    for n = 1:(metaDataSrt-1)
        [~] = fgetl(fid);
    end
    
    for n = 1 : metaDataEnd - metaDataSrt +1
        metaInfo{n,1} = fgetl(fid);
    end
else
    metaInfo = []; % Set meta info as null
end

% Redundant functions
% function timeVec = Time_String_To_Seconds(timeCell)
% %convert a cell array of time strings to the time in seconds
% N = numel(timeCell);
% timeVec = zeros(N,1);
% 
% tempTimeStr = regexp(timeCell{1},'\"', 'split');
% if length(tempTimeStr) == 1 % The time string does not contain "
%     timeCmptns = regexp(tempTimeStr{1},'\:', 'split');
%     if length(timeCmptns) == 2      % time is in M:S.S format
%         timeFormat = '%d:%d.%d';
%     elseif length(timeCmptns) == 3  % time is in H:M:S.S format
%         timeFormat = '%d:%d:%d.%d';
%     end
% else                        % The time string does contain "
%     timeCmptns = regexp(tempTimeStr{2},'\:', 'split');
%     if length(timeCmptns) == 2      % time is in M:S.S format
%         timeFormat = '="%d:%d.%d"';
%     elseif length(timeCmptns) == 3  % time is in H:M:S.S format
%         timeFormat = '="%d:%d:%d.%d"';
%     end
% end
% 
% for n = 1:N
%     timeElements = sscanf(timeCell{n},timeFormat);
%     if isempty(timeElements)
%         continue
%     end
%     % Last element of timeElements is the value after the decimal point.
%     % Convert back to decimal from units. The logic is if the last element
%     % is for example 55 this should be converted to 0.55. To do represent
%     % the number as 0.55E-2 and get the value of the exponent. Then divide
%     % the last element by 10^2.
%     expo = ceil(log10([timeElements(end)+ 0.1])); % The 0.1 is added to ensure that all 10^n values return n+1 as the exponent when rounding up
%     div = 10^expo;
%     if length(timeElements) == 3
%         multiplier = [60, 1, 1/div];
%     elseif length(timeElements) == 4
%         multiplier = [3600, 60, 1, 1/div];
%     end
%     if ~isempty(timeElements)
%         timeVec(n) = multiplier*timeElements;
%     end
% end
% 
% idx = find(diff(timeVec)<0);
% 
% if ~isempty(idx)
%     for ii = 1:length(idx)+1
%         if ii == 1
%             sIdx = 1;
%             eIdx = idx(ii);
%         elseif ii == length(idx)+1
%             sIdx = idx(ii-1)+1;
%             eIdx = N;
%         else
%             sIdx = idx(ii-1)+1;
%             eIdx = idx(ii);
%         end
%         timeVec(sIdx:eIdx) = timeVec(sIdx:eIdx)+3600*(ii-1);
%     end
% end

