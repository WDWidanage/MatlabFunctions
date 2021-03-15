function MaccorProfiles_Excel(fileName,matFileName,sheetNumber)
% Import Bitrode data and save the profiles in a mat file
%
% Input arguments 
%   fileName: String variable of Excel file to read
%
% Optional input arguments
%   matFileName: Name of mat file in which profiles are saved
%   sheetNumber: Integer value specifying the Excel sheet
%   
% W. D. Widanage 23/09/2013 (crunch crunch!)

try matFileName;
catch
    hourMinSec = fix(clock);
    matFileName = ['Bitrode_',date,'_',num2str(hourMinSec(4)),'_',num2str(hourMinSec(5)),'_',num2str(hourMinSec(6))];
end

try sheetNumber;
catch
    sheetNumber = 1;
end

data = xlsread(fileName,sheetNumber);

timeTmpVec = data(:,4)*3600;                % Time vector (s)
currentTmpVec = data(:,8);             % current (A) 
volTmpVec = data(:,9);                 % Voltage (V)
powTmpVec = data(:,10);                % Power (W)
capTmpVec = data(:,6);                % Capacity (Ah)
eneTmpVec = data(:,7);                % Energy (Wh)
temp1TmpVec = data(:,15);               % Temperature (degC)
temp2TmpVec = data(:,16);               % Temperature (degC)

% Save data in a mat file
save(matFileName,'timeTmpVec','currentTmpVec','volTmpVec','powTmpVec','capTmpVec','eneTmpVec','temp1TmpVec','temp2TmpVec')

end

