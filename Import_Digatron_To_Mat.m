function Import_Digatron_To_Mat(varargin)
% Syntax of use
%   Import_Bitrode_To_Mat()
%   Import_Bitrode_To_Mat(pwd)
%   Import_Bitrode_To_Mat(pwd,pwd)
%   Import_Bitrode_To_Mat(pwd,pwd,[headerRowPosition,dataRowPosition])
%   Import_Bitrode_To_Mat(pwd,pwd,[headerRowPosition,dataRowPosition],'saveNames',{matFileName1,...,matFileNameN})
%   Import_Bitrode_To_Mat('saveNames',{matFileName1,...,matFileNameN})
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
%   Input 4 (optional): A property-value input specifying file names to use
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

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
% validity - 09-July-2015
addOptional(parObj,'csvDetails',[], @checkCSVDetails);
addOptional(parObj,'savePath',[], @checkSavePath);
addOptional(parObj,'hdPos',[], @checkHdPos);
addParameter(parObj,'saveNames',[], @checkSaveNames);
addParameter(parObj,'metaPos',[], @checkMetaPos);


% Re-parse parObj - 09-July-2015
parse(parObj,varargin{:});


if isempty(parObj.Results.csvDetails)
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
    fullFileNames = Input_Process(parObj.Results.csvDetails);
end

%check all files are CSVs and exist.
fullFileNames = Check_Filenames(fullFileNames);

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
    data = Import_Data(currentFile,hdrData(1),hdrData(2));
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
function valid = checkCSVDetails(csvDetails)
% Check if csvDetails is a string or cell array of strings
if strcmp(csvDetails,'saveNames')
    valid = false;
elseif ischar(csvDetails)
    valid = true;
elseif all(cellfun(@ischar,csvDetails))
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




function fullFileNames = Input_Process(arg)
% parse first input to main function.
if ischar(arg)
    argCell = {arg}; %cell so it can be treated the same as if a cell array
elseif iscell(arg)
    argCell = arg;
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
        D = dir(sprintf('%s%s%s',filePath,filesep,'*.csv'));
        fileNames = {D.name}';
        if isempty(fileNames)
            error('No CSV files found in directory')
        else
            fullFileNames = cellfun(@(f) sprintf('%s%s%s',filePath,filesep,f),fileNames,'UniformOutput',false);
        end
        
    else
        %assume it is a file with path for now (check later)
        fullFileNames = argCell(1); %keep as cell array
    end
end





function fullFileNames = Check_Filenames(fullFileNames)
%check that the file names are CSVs and all exist.
for n = 1:numel(fullFileNames)
    [p,f,e] = fileparts(fullFileNames{n});
    
    %add CSV extension is there is none
    if isempty(e)
        e = '.csv';
    end
    
    %check extension is correct
    if ~strcmp(e,'.csv')
        error('%s is not a csv file.',fullFileNames{n});
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
    hdrDataVec = [30,35];
end






function dataStruct = Import_Data(fileName,rowColumnNames,rowDataStart)
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
    columnUnitText = fgetl(fid); % channel Units
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
[dataStruct,formatString] = Process_Channel_Names(columnNameText,columnUnitText,dataFirstLine);
nChannels = numel(dataStruct.channelNames);
dataStruct.csvInfo = dir(fileName);

%read in data (a 1 x nChannel cell array):
data = textscan(fid,formatString,'delimiter',',');
fclose('all');


% save columns of data to structure, depending on if string or numeric
dataNames = fieldnames(dataStruct.data);
for n = 1:nChannels
    if iscell(data{n})
        % Check whether the strings are time or not
        % If the cell string consisits of a ':' then this will be a time string column and therefore process the data 
        colonIdx = regexp(data{n}{1},':', 'once');
        %~isempty(strfind(lower(dataNames{n}),'time')) || ~isempty(strfind(lower(dataNames{n}),'channel_2'))
        if ~isempty(colonIdx)
%             dataStruct.data.(dataNames{n}) = Time_String_To_Seconds(data{n});  % Process time strings to seconds
        else
            dataStruct.data.(dataNames{n}) = char(data{n}); %char: requires much less memory than cell
        end
    else
        dataStruct.data.(dataNames{n}) = data{n};
        nanRows = isnan(data{n});
    end
end

% dataStruct.data = structfun(@(x) x(~nanRows,:),dataStruct.data,'UniformOutput',false); %remove any rows which has NaNs in the numeric columns
dataStruct.nSamples = size(dataStruct.data.(dataNames{n}),1);






function [dataStruct,formatString] = Process_Channel_Names(channelString,channelUnit,dataString)
%get channel names, and test to determine which channels are text and whcih
%are numeric.
%channelString (char array): row of spreadsheet containing channel names.
%dataString (char array): row of spreadsheet containing a row of data

% %use first line of data to see how many channels there are:
% nChannels = numel(strfind(dataString,','))+1;

%read all data channels as a string, then see which can be converted to
%numbers:
dataTextTmp = regexp(dataString,',','split');
empData = cellfun(@isempty,dataTextTmp);
dataText = dataTextTmp(~empData);
nChannels = numel(dataText);
dataConv = cellfun(@str2double,dataText);
formatString = repmat('%s',1,nChannels);
formatString(find(~isnan(dataConv))*2) = 'f'; %*2: every second character

%if there wasn't any column names row, then create a dummy one
if isempty(channelString)           % No headers are specified. Use default, Channel_1...Channel_N
    channelNamesTmp = sprintf('Channel_%d',1:nChannels);
    channelNames = regexp(channelNamesTmp,'Channel_[\d]{0,}','match'); % Get channel hdr names
    units = repmat({'None'},nChannels,1);
else
    
    expHdr = {'Step','Status','Step Time','Prog Time','Cycle','Cycle Level',...
              'Procedure','Voltage','Current','T_Bat_1','PresetU','PresetI',...
              'PMPresetI_1','PMPresetU_1','Capacity','AhCha','AhDch','AhStep',...
              'WhAccu','WhCha','WhDch','WhStep','Charge Factor','Power',...
              'Time Stamp','Watt','LogTemp001','LogTemp002','LogTemp0450104', 'AhAccu', 'AhPrev',...
              'PMMeasureI_1',	'PMMeasureU_1',	'V_DC_Link_1',	'UCapa_1',...
              'UTerminal_1','TEMP_MOS_1',	'12V_SUPPLY_1',	'5V_SUPPLY_1',	'TEMP_RACK_1','I_Raw_1'};    

    channelNamesTmp = regexp(channelString,'[\w\s-#]{3,}','match'); % Get channel hdr names
    
    count = 0;
    for ii = 1:length(channelNamesTmp)
        if ismember(channelNamesTmp{ii},expHdr)
            count = count + 1;
            channelNames{count} = channelNamesTmp{ii};
        end
    end
    

    % Handle units
    expHdrsWithUnits =  {'Voltage','Current','T_Bat_1','PresetU','PresetI',...
              'PMPresetI_1','PMPresetU_1','Capacity','AhCha','AhDch','AhStep',...
              'WhAccu','WhCha','WhDch','WhStep','ChargeFactor','Power',...
              'Time Stamp','Watt','LogTemp001','LogTemp002','LogTemp0450104','AhAccu',}; 
    nChannels = length(channelNames);
    channelText_Units_Tmp = split(channelUnit,',');
    unitsTmp = matlab.lang.makeValidName(channelText_Units_Tmp);
    units = unitsTmp(1:end-1)';
    
end

% channelText = regexp(channelString,'(?<=")[\w\s,#-:]{2,}(?=")','match');

%create a structure to store to the data in:
structNames = matlab.lang.makeValidName(channelNames); %make channel names valid field names
structDummy = [structNames;cell(1,nChannels)];
dataStruct.data = struct(structDummy{:});

%add channel info to structure
dataStruct.channelNames = channelNames;
dataStruct.channelUnits = units;





function timeVec = Time_String_To_Seconds(timeCell)
%convert a cell array of time strings to the time in seconds
N = numel(timeCell);
timeVec = zeros(N,1);

tempTimeStr = regexp(timeCell{1},'\"', 'split');
if length(tempTimeStr) == 1 % The time string does not contain "
    timeCmptns = regexp(tempTimeStr{1},'\:', 'split');
    if length(timeCmptns) == 2      % time is in M:S.S format
        timeFormat = '%d:%d.%d';
    elseif length(timeCmptns) == 3  % time is in H:M:S.S format
        timeFormat = '%d:%d:%d.%d';
    end
else                        % The time string does contain "
    timeCmptns = regexp(tempTimeStr{2},'\:', 'split');
    if length(timeCmptns) == 2      % time is in M:S.S format
        timeFormat = '="%d:%d.%d"';
    elseif length(timeCmptns) == 3  % time is in H:M:S.S format
        timeFormat = '="%d:%d:%d.%d"';
    end
end

for n = 1:N
    timeElements = sscanf(timeCell{n},timeFormat);
    if isempty(timeElements)
        continue
    end
    % Last element of timeElements is the value after the decimal point.
    % Convert back to decimal from units. The logic is if the last element
    % is for example 55 this should be converted to 0.55. To do represent
    % the number as 0.55E-2 and get the value of the exponent. Then divide
    % the last element by 10^2.
    expo = ceil(log10([timeElements(end)+ 0.1])); % The 0.1 is added to ensure that all 10^n values return n+1 as the exponent when rounding up
    div = 10^expo;
    if length(timeElements) == 3
        multiplier = [60, 1, 1/div];
    elseif length(timeElements) == 4
        multiplier = [3600, 60, 1, 1/div];
    end
    if ~isempty(timeElements)
        timeVec(n) = multiplier*timeElements;    
    end
end

idx = find(diff(timeVec)<0);

if ~isempty(idx)
    for ii = 1:length(idx)+1
        if ii == 1
            sIdx = 1;
            eIdx = idx(ii);
        elseif ii == length(idx)+1
            sIdx = idx(ii-1)+1;
            eIdx = N;
        else
            sIdx = idx(ii-1)+1;
            eIdx = idx(ii);
        end
        timeVec(sIdx:eIdx) = timeVec(sIdx:eIdx)+3600*(ii-1);
    end
end