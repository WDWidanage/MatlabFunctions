% Sanity check of data
%
% W.D. Widanage (15/10/2018)
function maccorDataPlot(data)

try
    timeAxis = data.TestTime_s_/3600; % Time [hours]
catch
    timeAxis = data.TestTime/3600; % Time [hours]
end
try
    currSig = data.Amps; % Cuurrent [A]
catch
    currSig = data.Current_A_;
end
try
    volSig = data.Volts; % Voltage [V]
catch
    volSig = data.Voltage_V_;
end
try
    tempSig = data.Temp1; % Temperature [degC]
catch
    tempSig = data.VAR1; %data.AUX1_C_;
end
try
    capSig = data.Amp_hr; % Capacity [Ah]
catch
    capSig = data.Cap__Ah_;
end

stepNums = data.Step;
timeAxis = timeAxis - timeAxis(1);

figure
ax1 = subplot(6,1,1);
plot(timeAxis,currSig,'. -');
xlabel('Time (Hours)'); ylabel('Current (A)');
ax2 = subplot(6,1,2);
plot(timeAxis,volSig,'. -');
xlabel('Time (Hours)'); ylabel('Voltage (V)');
ax3 = subplot(6,1,3);
plot(timeAxis,volSig.*currSig,'. -');
xlabel('Time (Hours)'); ylabel('Power (W)');
ax4 = subplot(6,1,4);
plot(timeAxis,capSig,'. -');
xlabel('Time (Hours)'); ylabel({'Charge','throughput (Ah)'});
ax5 = subplot(6,1,5);
plot(timeAxis,tempSig,'. -');
xlabel('Time (Hours)'); ylabel({'Temperature','(degC)'});
ax6 = subplot(6,1,6);
plot(timeAxis,stepNums,'. -');
xlabel('Time (Hours)'); ylabel({'Step number','(-)'});


linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
end

