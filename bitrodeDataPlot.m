% Sanity check of data
%
% W.D. Widanage (05/03/2019)
function bitrodeDataPlot(data)

timeAxis = data.TotalTime/3600; % Time [hours]

currSig = data.Current; % Cuurrent [A]
volSig = data.Voltage; % Voltage [V]
try
    tempSig = data.TemperatureA1; % Temperature [degC]
catch
    tempSig = zeros(size(data.Current)); % Temperature [degC]
end
capSig = data.Amp_Hours; % Capacity [Ah]

figure
ax1 = subplot(5,1,1);
plot(timeAxis,currSig,'. -');
xlabel('Time (Hours)'); ylabel('Current (A)');
ax2 = subplot(5,1,2);
plot(timeAxis,volSig,'. -');
xlabel('Time (Hours)'); ylabel('Voltage (V)');
ax3 = subplot(5,1,3);
plot(timeAxis,volSig.*currSig,'. -');
xlabel('Time (Hours)'); ylabel('Power (W)');
ax4 = subplot(5,1,4);
plot(timeAxis,capSig,'. -');
xlabel('Time (Hours)'); ylabel({'Charge','throughput (Ah)'});
ax5 = subplot(5,1,5);
plot(timeAxis,tempSig,'. -');
xlabel('Time (Hours)'); ylabel({'Temperature','(degC)'});


linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
end