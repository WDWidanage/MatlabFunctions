function [ssResults] = snapShotAnalysis_Sam3Ah(data,varargin)
% Extract capacity and resitances of the Samsung 3Ah cell.
% The data must be collected from a Maccor tester
%
% W.D. Widanage 04/11/18 (Switch)

p = inputParser;

% Create variable names and assign default values after checking the value
addRequired(p,'data');

% Optional parameteres. Indices of capacity and resistances
addParameter(p,'posCapC25',[5,2]);
addParameter(p,'posCap1C',[3,2]);
addParameter(p,'posRes2C',[7,1]); % 31s pulse
addParameter(p,'posRes3C',[9,1]); % 10s pulse

% Re-parse parObj
parse(p,data,varargin{:});

idxExc = find_exc_segments(data.Current,'timeVec',data.TotalTime,'plotSeg',1);

% C/3 capacity
idxTmp0 = p.Results.posCapC25(1); idxTmp1 = p.Results.posCapC25(2);
ssResults.Capacity.C3 = data.Amp_hr(idxExc(idxTmp0,idxTmp1));

% 1C capacity
idxTmp0 = p.Results.posCap1C(1); idxTmp1 = p.Results.posCap1C(2);
ssResults.Capacity.C1 = data.Amp_hr(idxExc(idxTmp0,idxTmp1));

% 2C pulse resistances
idxTmp0 = p.Results.posRes2C(1); idxTmp1 = p.Results.posRes2C(2);
startIdx = idxExc(idxTmp0,idxTmp1); endIdx = idxExc(idxTmp0,2);
R0 = (data.Volts(startIdx+1) - data.Volts(startIdx))/data.Amps(startIdx+1);
R35 = (data.Volts(endIdx) - data.Volts(startIdx))/data.Amps(startIdx+1);
ssResults.Resistance.Pulse_2C.R0 = R0;
ssResults.Resistance.Pulse_2C.R35 = R35;


% 3C pulse resistances
idxTmp0 = p.Results.posRes3C(1); idxTmp1 = p.Results.posRes3C(2);
startIdx = idxExc(idxTmp0,idxTmp1); endIdx = idxExc(idxTmp0,2);
R0 = (data.Volts(startIdx+1) - data.Volts(startIdx))/data.Amps(startIdx+1);
R10 = (data.Volts(endIdx) - data.Volts(startIdx))/data.Amps(startIdx+1);
ssResults.Resistance.Pulse_3C.R0 = R0;
ssResults.Resistance.Pulse_3C.R10 = R10;

% Temperatures before discharge or pulse
% C3 capacity
idxTmp0 = p.Results.posCapC25(1);
ssResults.Temperature.C3 = data.Temp1(idxExc(idxTmp0,1));

% 1C capacity
idxTmp0 = p.Results.posCap1C(1);
ssResults.Temperature.C1 = data.Temp1(idxExc(idxTmp0,1));

% 2C resistance
idxTmp0 = p.Results.posRes2C(1); idxTmp1 = p.Results.posRes2C(2);
startIdx = idxExc(idxTmp0,idxTmp1); 
ssResults.Temperature.Pulse_2C = data.Temp1(startIdx);

% 3C resistance
idxTmp0 = p.Results.posRes3C(1); idxTmp1 = p.Results.posRes3C(2);
startIdx = idxExc(idxTmp0,idxTmp1); 
ssResults.Temperature.Pulse_3C = data.Temp1(startIdx);

end

