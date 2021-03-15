function ImportProfiles_Excel(fileName,matFileName,sheetNumber,convertTime)
% Import Bitrode data and save the profiles in a mat file
%
% Input arguments 
%   fileName: String variable of Excel file to read
%
% Optional input arguments
%   matFileName: Name of mat file in which profiles are saved
%   sheetNumber: Integer value specifying the Excel sheet
%   convertTime: Interger value of 1 if time vector is in H:M:S.S format
%   from Bitrode or 0 otherwise Default convertTime = 0;
%   
% W. D. Widanage 23/09/2013 (crunch crunch!)
% Edited to save data as a structure variable 09/07/2015

try matFileName;
catch
    hourMinSec = fix(clock);
    matFileName = ['Bitrode_',date,'_',num2str(hourMinSec(4)),'_',num2str(hourMinSec(5)),'_',num2str(hourMinSec(6))];
end

try sheetNumber;
catch
    sheetNumber = 1;
end

try convertTime;
    if convertTime == 1
        multiplier = 186400;
    else 
        multiplier = 1;
    end
catch
    multiplier = 1;
end

d = xlsread(fileName,sheetNumber);
[nSamples, nColumns] = size(d);


timeVec = multiplier*d(:,1)-multiplier*d(1,1);  % Total time
Cycle = d(:,2);                % Cycle
LoopCounter_1 = d(:,3);        % Loop counter 1
LoopCounter_2 = d(:,4);        % Loop counter 2
LoopCounter_3 = d(:,5);        % Loop counter 3
Step = d(:,6);                 % Step
StepTime = multiplier*d(:,7);             % Step Time
currentVec = d(:,8);               % current (A) 
volVec = d(:,9);               % Voltage (V)
powVec = d(:,10);                % Power (W)
capVec = d(:,11);            % Capacity (Ah)
eneVec = d(:,12);           % Energy (Wh)
tempVec = d(:,13);               % Temperature (degC)

% Check if filename is passed with .xlsm extension
ext = regexp(fileName,'.xlsm','match');
if isempty(ext)
    fileName = [fileName,'.xlsm'];
end

channelUnits = {'S';'None';'None';'None';'None';'None';'S';'A';'V';'W';'AH';'WH';'degC'};
excelInfo = dir(fileName);
% Save data in a mat file
save(matFileName,'timeVec','currentVec','volVec','capVec','tempVec','eneVec','powVec');

end

