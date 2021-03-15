% Sanity check of data
%
% W.D. Widanage (17/07/2019)
function vmp3DataPlot(data)

timeAxis = data.time_s/3600; % Time [hours]

currSig = data.ave_I_mA; % Cuurrent [A]
volSig = data.ave_Ewe_V; % Voltage [V]

freqSig = data.freq_Hz;     % Frequency [Hz]
magZ = data.mag_Z_Ohm;      % Impedance magnitude [Ohms]
phaseZ = data.Phase_Z_deg;  % Impedance phase [deg]
reZ = data.Re_Z_Ohm;        % Impedance real [Ohms]
imZ = data.NAME;            % Impedance imag [Ohms]

figure
ax1 = subplot(2,1,1);
plot(timeAxis,currSig,'. -');
xlabel('Time (Hours)'); ylabel('Current (A)');
ax2 = subplot(2,1,2);
plot(timeAxis,volSig,'. -');
xlabel('Time (Hours)'); ylabel('Voltage (V)');

figure
ax3 = subplot(2,1,1);
plot(freqSig,magZ,'. -');
xlabel('Frequency (Hz)'); ylabel('Magnitude (Ohms)');
ax4 = subplot(2,1,2);
plot(freqSig,phaseZ,'. -');
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');

figure
plot(reZ,imZ,'. -');
xlabel('Real Z(Ohms)'); ylabel('-Imaginary Z (Ohms)');



linkaxes([ax1,ax2],'x')
linkaxes([ax3,ax4],'x')

end