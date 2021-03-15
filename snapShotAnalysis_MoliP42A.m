function ssResults = snapShotAnalysis_MoliP42A(data,varargin)
% Extract capacity, resitances of the Moli P42A cell.
% The data must be collected from a Digatron tester and have followed the
% Moli 42A LTA Reference Preformance Test protocol.
%
% W.D. Widanage 06/12/20 (Quietness)

% Analysis results
%  Capacity
%       C/3
%       Residual
%
% Resistances
%      100% C3, 2C  30s
%       95% C3, 2C  30s
%       90% C3, 2C  30s
%       80% C3, 2C  30s
%       70% C3, 2C  30s
%       50% C3, 2C  30s
%       30% C3, 2C, 30s
%       20% C3, 2C, 30S
%       10% C3, 2C, 30s
%        5% C3, 2C, 30s

p = inputParser;

% Create variable names and assign default values after checking the value
addRequired(p,'data');

% Optional parameteres. Indices of capacity and resistances
addParameter(p,'capC3Step',13);
addParameter(p,'capC3ResStep',16);

addParameter(p,'resC3DSteps',2); % Step numbers for C/3 discharge pulses
addParameter(p,'res2CDSteps',8); % Step numbers for 2C discharge pulses
addParameter(p,'resC3CSteps',5);             % Step number for C/3 charge pulse
addParameter(p,'res2CCSteps',11); % Step numbers for 2C charge pulses
addParameter(p,'resTimes',[1,10]);              % Resistance time

% Re-parse parObj
parse(p,data,varargin{:});

capC3Step = p.Results.capC3Step;
capC3ResStep = p.Results.capC3ResStep;
resC3DSteps = p.Results.resC3DSteps;
res2CDSteps = p.Results.res2CDSteps;
resC3CSteps = p.Results.resC3CSteps;
res2CCSteps = p.Results.res2CCSteps;
resTimes = p.Results.resTimes;

Cn = 4.2; % Nominal capacity [Ah]

% Capacity
dataTmpC3 = extractStepSegments(data,data.Step,capC3Step);
[~,aHC3] = CoulombCounting(abs(dataTmpC3(1).Current),dataTmpC3(1).ProgTime/3600,Cn,Cn,0);

dataTmpRes = extractStepSegments(data,data.Step,capC3ResStep);
[~,aHC3Res] = CoulombCounting(abs(dataTmpRes(1).Current),dataTmpRes(1).ProgTime/3600,Cn,Cn,0);

capC3 = aHC3(end);
capRes = aHC3Res(end);
ssResults.Capacity.C3_Ah = capC3;
ssResults.Capacity.C3res_Ah = capC3 + capRes;


% Resistances: C/3 discharges
resSoC_D = [100,95,90,80,70,50,30,20,10,5];       % Refrence soc points for discharge pulses
dataC3DPlTmp1 = extractStepSegments(data,data.Step,20);
dataC3DPlTmp2 = extractStepSegments(data,data.Step,resC3DSteps);
dataC3DPulse = [dataC3DPlTmp1;dataC3DPlTmp2];
for ii = 1:length(dataC3DPulse)
    if isnan(resC3DSteps)
            ssResults.Resistance.Discharge.(sprintf('resC3_30s_%dper_Ohms',resSoC_D(ii))) = NaN;
            ssResults.Resistance.Discharge.soc(ii,1) = resSoC_D(ii);
            ssResults.Resistance.Discharge.resC3_30s_Ohms(ii,1) = NaN;
            for rr = 1:length(resTimes)
                ssResults.Resistance.Discharge.Pulse2C.(['SoC',num2str(resSoC_D(ii))]).(['R',resStr]) = NaN;
            end
    else 
    dataTmpVol = dataC3DPulse(ii);
    [resTmp,opt] = dcirDigatron(dataTmpVol,'resTimes',resTimes,'Cn',Cn);
    resTmp = resTmp*1000;
    for rr = 1:length(resTimes)
        if resTimes(rr)<1
            resStr = ['0p',num2str(10*resTimes(rr))];
        else
            resStr = num2str(resTimes(rr));
        end
        ssResults.Resistance.Discharge.PulseC3.(['SoC',num2str(resSoC_D(ii))]).(['R',resStr,'mOhms']) = resTmp(rr);
        ssResults.Resistance.Discharge.PulseC3.(['SoC',num2str(resSoC_D(ii))]).(['R',resStr,'mOhms_Msg']) = opt.msg;
    end
    
    deltaV_30 = dataTmpVol.Voltage(1) - dataTmpVol.Voltage(end);
    appI = abs(dataTmpVol.Current(end));
    ssResults.Resistance.Discharge.(sprintf('resC3_30s_%dper_Ohms',resSoC_D(ii))) = deltaV_30/appI;
    ssResults.Resistance.Discharge.soc(ii,1) = resSoC_D(ii);
    ssResults.Resistance.Discharge.resC3_30s_Ohms(ii,1) = deltaV_30/appI;
    dataTmpVol = [];
    end
end

% Resistances: 2C discharges
resSoC_D = [95,90,80,70,50,30,20,10,5];       % Refrence soc points for discharge pulses
data2CDPlTmp = extractStepSegments(data,data.Step,res2CDSteps);
data2CDPulse = data2CDPlTmp(2:end);
for ii = 1:length(data2CDPulse)
    if isnan(res2CDSteps)
            ssResults.Resistance.Discharge.(sprintf('res2C_30s_%dper_Ohms',resSoC_D(ii))) = NaN;
            ssResults.Resistance.Discharge.soc(ii,1) = resSoC_D(ii);
            ssResults.Resistance.Discharge.res2C_10s_Ohms(ii,1) = NaN;
            for rr = 1:length(resTimes)
                ssResults.Resistance.Discharge.Pulse2C.(['SoC',num2str(resSoC_D(ii))]).(['R',resStr]) = NaN;
            end
    else 
    dataTmpVol = data2CDPulse(ii);
    [resTmp,opt] = dcirDigatron(dataTmpVol,'resTimes',resTimes,'Cn',Cn);
    resTmp = resTmp*1000;
    for rr = 1:length(resTimes)
        if resTimes(rr)<1
            resStr = ['0p',num2str(10*resTimes(rr))];
        else
            resStr = num2str(resTimes(rr));
        end
        ssResults.Resistance.Discharge.Pulse2C.(['SoC',num2str(resSoC_D(ii))]).(['R',resStr,'mOhms']) = resTmp(rr);
        ssResults.Resistance.Discharge.Pulse2C.(['SoC',num2str(resSoC_D(ii))]).(['R',resStr,'mOhms_Msg']) = opt.msg;
    end
    
    deltaV_30 = dataTmpVol.Voltage(1) - dataTmpVol.Voltage(end);
    appI = abs(dataTmpVol.Current(end));
    ssResults.Resistance.Discharge.(sprintf('res2C_30s_%dper_Ohms',resSoC_D(ii))) = deltaV_30/appI;
    ssResults.Resistance.Discharge.soc(ii,1) = resSoC_D(ii);
    ssResults.Resistance.Discharge.res2C_30s_Ohms(ii,1) = deltaV_30/appI;
    dataTmpVol = [];
    end
end

% Resistance: C/3 charge
resSoC_C = [95,90,80,70,50,30,20,10,5];          % Refrence soc points for charge pulses
dataC3CPulseTmp = extractStepSegments(data,data.Step,resC3CSteps);
if length(dataC3CPulseTmp) > 9
    dataC3CPulse = dataC3CPulseTmp(2:end);
else
    dataC3CPulse = dataC3CPulseTmp;
end
for ii = 1:length(dataC3CPulse)
    if isnan(resC3CSteps)
            ssResults.Resistance.Charge.(sprintf('resC3_30s_%dper_Ohms',resSoC_C(ii))) = NaN;
            ssResults.Resistance.Charge.soc(ii,1) = resSoC_C(ii);
            ssResults.Resistance.Charge.resC3_30s_Ohms(ii,1) = NaN;
            for rr = 1:length(resTimes)
                ssResults.Resistance.Charge.PulseC3.(['SoC',num2str(resSoC_C(ii))]).(['R',resStr]) = NaN;
            end
    else 
    dataTmpVol = dataC3CPulse(ii);
    [resTmp,opt] = dcirDigatron(dataTmpVol,'resTimes',resTimes,'Cn',Cn);
    resTmp = resTmp*1000;
    for rr = 1:length(resTimes)
        if resTimes(rr)<1
            resStr = ['0p',num2str(10*resTimes(rr))];
        else
            resStr = num2str(resTimes(rr));
        end
        ssResults.Resistance.Charge.PulseC3.(['SoC',num2str(resSoC_C(ii))]).(['R',resStr,'mOhms']) = resTmp(rr);
        ssResults.Resistance.Charge.PulseC3.(['SoC',num2str(resSoC_C(ii))]).(['R',resStr,'mOhms_Msg']) = opt.msg;
    end
    
    deltaV_30 = abs(dataTmpVol.Voltage(1) - dataTmpVol.Voltage(end));
    appI = abs(dataTmpVol.Current(end));
    ssResults.Resistance.Charge.(sprintf('resC3_30s_%dper_Ohms',resSoC_C(ii))) = deltaV_30/appI;
    ssResults.Resistance.Charge.soc(ii,1) = resSoC_C(ii);
    ssResults.Resistance.Charge.resC3_30s_Ohms(ii,1) = deltaV_30/appI;
    dataTmpVol = [];
    end
end

% Resistance: 2C charge
resSoC_C = [95,90,80,70,50,30,20,10,5];          % Refrence soc points for charge pulses
data2CCPlTmp = extractStepSegments(data,data.Step,res2CCSteps);
if length(data2CCPlTmp)>9
    data2CCPulse = data2CCPlTmp(2:end);
else
    data2CCPulse = data2CCPlTmp;
end
for ii = 1:length(data2CCPulse)
    if isnan(res2CCSteps)
            ssResults.Resistance.Charge.(sprintf('res2C_30s_%dper_Ohms',resSoC_C(ii))) = NaN;
            ssResults.Resistance.Charge.soc(ii,1) = resSoC_C(ii);
            ssResults.Resistance.Charge.res2C_30s_Ohms(ii,1) = NaN;
            for rr = 1:length(resTimes)
                ssResults.Resistance.Charge.Pulse2C.(['SoC',num2str(resSoC_C(ii))]).(['R',resStr]) = NaN;
            end
    else 
    dataTmpVol = data2CCPulse(ii);
    [resTmp,opt] = dcirDigatron(dataTmpVol,'resTimes',resTimes,'Cn',Cn);
    resTmp = resTmp*1000;
    for rr = 1:length(resTimes)
        if resTimes(rr)<1
            resStr = ['0p',num2str(10*resTimes(rr))];
        else
            resStr = num2str(resTimes(rr));
        end
        ssResults.Resistance.Charge.Pulse2C.(['SoC',num2str(resSoC_C(ii))]).(['R',resStr,'mOhms']) = resTmp(rr);
        ssResults.Resistance.Charge.Pulse2C.(['SoC',num2str(resSoC_C(ii))]).(['R',resStr,'mOhms_Msg']) = opt.msg;
    end
    
    deltaV_30 = abs(dataTmpVol.Voltage(1) - dataTmpVol.Voltage(end));
    appI = abs(dataTmpVol.Current(end));
    ssResults.Resistance.Charge.(sprintf('res2C_30s_%dper_Ohms',resSoC_C(ii))) = deltaV_30/appI;
    ssResults.Resistance.Charge.soc(ii,1) = resSoC_C(ii);
    ssResults.Resistance.Charge.res2C_30s_Ohms(ii,1) = deltaV_30/appI;
    dataTmpVol = [];
    end
end
end
