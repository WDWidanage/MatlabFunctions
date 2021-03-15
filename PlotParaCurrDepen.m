function  PlotParaCurrDepen(cratesCh,cratesDh,socHPPC,tempBP,paraResutlsCh,paraResutlsDh, cellNo,figOffSet)
%
% Plot the ECM parameter variations as a function of current for all SoCs
% and temperatures
%
% Input arguments:
%   cratesCh: Charge current matrix          size: number of pulses x SoC values x temperature
%   cratesDh: Discharge current matrix       size: number of pulses x SoC values x temperature
%   socHPPC: SoC break points matrix           size: SoC points x cell numbers x temperature
%   tempBP: Temperature break points         size: Temp points x 1
%   paraResutlsCh: Charge parameters of ECM cell variable  size: Number of cells x charge pulses x SoC values x temperature
%   paraResutlsDh: Discharge parameters of ECM cell variable  size: Number of cells x charge pulses x SoC values x temperature
%   cellNo: Indicate which of the cell results to be plots. Scalar value
%           and should not exceed first dimension of parameter matrices
%
% W. D. Widanage 30/04/2014 (Crunch crunch!)


[numCells,numPulses, numSoCPts, numTemp] = size(paraResutlsCh);

if cellNo > numCells
    error('fcn:cellArgChk','Specified cell exceeds number of cells over which parameters are estimated')
end

try figOffSet
catch
    figOffSet = 0;
end

% Extract parameter set for specified cell
tmpParaC = squeeze(paraResutlsCh(cellNo,:,:,:));    % Cell variable of size pulses x soc values x temperature values
tmpParaD = squeeze(paraResutlsDh(cellNo,:,:,:));    % Cell variable of size pulses x soc values x temperature values

socTemp = squeeze(socHPPC(:,cellNo,:));             % Extract SoC for given cell decreasing SoC x Temperature
SoCBP = mean(socTemp,2);

for tt = 1:numTemp
    for zz = 1:numSoCPts
        
        % Get the corresponding temperature and SoC
        temperature = tempBP(tt);
        SoC = SoCBP(zz);
        
        
        chargeCurrents = cratesCh(:,numSoCPts-zz+1,tt);
        dischargeCurrents = cratesDh(:,numSoCPts-zz+1,tt);
        
        
        for pp = 1:numPulses
            
            % Charge parameters
            paraCPulses = tmpParaC(:,numSoCPts-zz+1,tt);  % Cell variable contains all the charge parameters for a given cell, all pulses, given temperature and given soc. cell size:  numPulses x 1
            paraC(:,pp) = paraCPulses{pp};
            
            % Disharge parameters
            paraDPulses = tmpParaD(:,numSoCPts-zz+1,tt);  % Cell variable contains all the discharge parameters for a given cell, all pulses, given temperature and given soc. cell size: numPulses x 1
            paraD(:,pp) = paraDPulses{pp};
            
        end
        % Get order of ECM
        [orderTmp,~] = size(paraC);
        order = (orderTmp-2)/2;
        
        % Flip and join current to a single varaible
        currents = [flipud(chargeCurrents);dischargeCurrents];
        para = [fliplr(paraC),paraD];
               
        % Plots are done for the Ro, TF gain G = Ro+Rp and time constants
        Ro = para(1,:);
        Rp = para(3:3+order-1,:);
        G = sum([Ro;Rp]);
        tau = para(3+order:end,:);
        
        figure(figOffSet); hold on;
        subplot(numTemp,1,tt)
        plot(currents,Ro,'. -');
        ylabel('Ro (Ohms)'); title([num2str(tempBP(tt)),'degC']); legend(['SoC ',num2str(SoC),'%']);
        
        figure(figOffSet+1); hold on;
        subplot(numTemp,1,tt)
        plot(currents,G,'. -');
        ylabel('Gain (Ohms)'); title([num2str(tempBP(tt)),'degC'])
        
        for par = 1:order            
            figure(figOffSet+par+1); hold on;
            subplot(numTemp,1,tt)
            plot(currents,tau(par,:),'. -');
            ylabel('Tau (s)'); title([num2str(tempBP(tt)),'degC'])
        end
        
    end
end
end

