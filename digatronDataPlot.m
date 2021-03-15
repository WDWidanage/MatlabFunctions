function digatronDataPlot(data,varargin)
% Sanity check of data
%
% W.D. Widanage (15/10/2018)

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments

addParameter(parObj,'resetTime','n')
addParameter(parObj,'timeUnit','H')

parse(parObj,varargin{:});

if ismember(parObj.Results.resetTime,{'y','Y'})
    timeAxis = (data.ProgTime - data.ProgTime(1)); % Time [s]
    xlabelStr = 'Time (Seconds)';
else
    timeAxis = data.ProgTime; % Time [s]
    xlabelStr = 'Time (Seconds)';
end

if ismember(parObj.Results.timeUnit,{'H'})
    timeAxis = timeAxis/3600;
    xlabelStr = 'Time (Hours)';
end

currSig = data.Current; % Cuurrent [A]
volSig = data.Voltage; % Voltage [V]
try
    tempSig = data.LogTemp001; % Temperature [degC]
catch
    tempSig = [data.LogTempPositive,data.LogTempMid,data.LogTempNegative]; % Temperature [degC]
end

capSig = data.AhAccu; % Capacity [Ah]
% stepNum = data.Step(1:end-2);

figure
ax1 = subplot(6,1,1);
plot(timeAxis,currSig,'. -');
xlabel(xlabelStr); ylabel('Current (A)');
ax2 = subplot(6,1,2);
plot(timeAxis,volSig,'. -');
xlabel(xlabelStr); ylabel('Voltage (V)');
ax3 = subplot(6,1,3);
plot(timeAxis,data.WhAccu,'. -');
xlabel(xlabelStr); ylabel('WhAccu (Wh)');
ax4 = subplot(6,1,4);
plot(timeAxis,capSig,'. -');
xlabel(xlabelStr); ylabel({'Charge','throughput (Ah)'});
ax5 = subplot(6,1,5);
plot(timeAxis,tempSig,'. -');
% ylim([22,30])
xlabel(xlabelStr); ylabel({'Temperature','(degC)'});
ax6 = subplot(6,1,6);
plot(timeAxis(1:end-5),data.Step(1:end-5),'. -');
xlabel(xlabelStr); ylabel({'Step num','(-)'});

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
end