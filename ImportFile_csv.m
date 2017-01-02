function ImportFile_csv(filename, Nheader)
% Imports data from the specified file, used initially for importing csv
% files from PEC battery tester for Project SuperLIB

% Input arguments
%  filename:  string variable of file to read, 'myfile.csv' or 'C:\my.csv'
%  Nheader: an interger value specifiyong how many header lines to leave
%           out before reading the data
%
% Usage
%  ImportFile_csv('myfile.csv');
%
% Creates a matrix varaible called data in the workspace with the
% extracted numeric values.
%
% Author
%  W.D. Widanage 05-08-2012 (Progress!!)

Delimiter = ',';

% Import the file
newData = importdata(filename, Delimiter, Nheader);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData);
for i = 1:length(vars)
    assignin('base', vars{i}, newData.(vars{i}));
end

