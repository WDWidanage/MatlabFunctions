function [currcd,volcd] = currentdependency(tempValues,tempBP,cellNumbers,hPPCStartMatrix,startEndRegions,Cn)
% Function to calculate current dependency
%
% Input arguments from data_analye_HPPC
%   tempValues       - Temperatures of interest
%   tempBP           - Temperatures at which experiments are performed
%   cellNumbers      - Cell identifier
%   hPPCStartMatrix  - Temperature against cell number
%   startEndRegions  - Define the main discharge pulses and end of HPPC pulse. These pulse excitation regions are used to determing total capacity extracted over the HPPC test
%   Cn               - To be used as the normalising value based on capacity extracts from HPPC test: see varaible cumCap once program is ru
%
% Output arguments:
%   currcd           - mean current for each pulse
%   voldc            - corresponding max voltage response
%   function plots currcd vs volcd
%
% T. R. B. Grandjean 29/04/2014

%%
volcd  = NaN([length(tempValues)*4*length(cellNumbers),10]) ; % initialise voltage array for plots, 4 SOCS, 10 pulses
currcd = NaN([length(tempValues)*4*length(cellNumbers),10]) ; % initialise current array for plots, 4 SOCS, 10 pulses
for ii = 1: length(tempValues) % loop for each temperature 

    [dummy,idxTemperature] = ismember(tempValues(ii),tempBP);
   
    for jj = 1:length(cellNumbers) % loop for cells
     
        % load data
        if tempValues(ii)<0
            tmp = ['m',num2str(abs(tempValues(ii)))];
            fileLocation = [tmp,'degC'];
        else
            fileLocation = [num2str(tempValues(ii)),'degC'];
        end
        matFileName = ['Profile_HPPC_Cell',num2str(cellNumbers(jj)),'_',num2str(tempValues(ii)),'degC'];
        load([fileLocation,'\',matFileName])
        
        timeHours = timeVec/3600;   % Convert time vector into hours
        currVec = -currentVec;      % -ve is charge
        
        % Obtain indices of the excitation and relaxation test portion
        [~,idx_ExcRelax] = find_exc_segments(currVec,1);
        HPPCStart = hPPCStartMatrix{idxTemperature,cellNumbers(jj)};
        
        if ~isempty(HPPCStart)
            % Obtain cumulative capacity from Bitrode tester to calculate SoC.
            % Capacity is calulcated from when fully charged to fully discharged
            for cc = 1:length(startEndRegions)
                pp = startEndRegions(cc);
                idxTmp = idx_ExcRelax(pp,1):idx_ExcRelax(pp,2);
                capTmp = capVec(idxTmp);            % Extract corresponding capacity signal
                maxCapDel = max(abs(capTmp));       % Get maximum change in Ah to calculate change in SoC
                % Store cummulative charge and discharge capacity for each cell
                % at each temperature
                if cc == 1
                    cumCap(cc,jj,ii) = maxCapDel;
                else
                    cumCap(cc,jj,ii) = cumCap(cc-1,jj,ii) + maxCapDel;
                end
            end
            normCap = Cn; % max(cumCap(:,jj,ii));
            socEstimate =  (normCap - cumCap(:,jj,ii))./normCap;
            [~,idx_soc_HPPCStart]= ismember(HPPCStart - 1,startEndRegions);    % Get the indices prior to HPPC cycle start to extract SoC at that point
            socHPPC(:,jj,ii) = socEstimate(idx_soc_HPPCStart)*100;            % The correponding soc prior to each HPPC cycle for each cell and temperature
            
            for kk = 1:length(HPPCStart) % loop for each SOC 
                
                % Check number of Pulses = 10
                if length(idx_ExcRelax)+1 - HPPCStart(kk) < 12
                    nPulses = floor((length(idx_ExcRelax)+1 - HPPCStart(kk))/2);
                else
                    nPulses = 5;
                end
               
                % Select HPPC pulses for a given SoC
                ps = HPPCStart(kk);
                idxTmp = idx_ExcRelax(ps,1):idx_ExcRelax(ps+2*nPulses-1,2);
                timeTmp = timeVec(idxTmp);
                currTmp = currVec(idxTmp);
                volTmp = volVec(idxTmp);

                    % Obtain indices of the excitation and relaxation test portion
                    [~,idx_Pulse] = find_exc_segments(currTmp,1);             
                    for p = 1:nPulses  % loop for each pulse
                        
                        cidx = 2*(nPulses+1-p) ; % charge pulse select: 10, 8, 6, 4, 2
                        didx = 2*p-1 ;           % discharge pulse select: 1, 3, 5, 7, 9
                        
                        % check number of sample points on current
                        % for charge pulses
                        cdatapts = 1 ; 
                        while abs(currTmp(idx_Pulse(cidx)+cdatapts+1) - currTmp(idx_Pulse(cidx))) > 1
                            cdatapts = cdatapts + 1 ;
                        end
                        
                        % for discharge pulses
                        ddatapts = 1 ; 
                        while abs(currTmp(idx_Pulse(didx)+ddatapts+1) - currTmp(idx_Pulse(didx))) > 1 
                            ddatapts = ddatapts + 1 ;
                        end

                        % compute max voltage and average current
                        datacount = ((ii-1)*length(cellNumbers)+(jj-1))*length(HPPCStart)+kk; % row index for voltage and current arrays, (SOC, cell, temp) x (pulses)
                        % charge pulses
                        volcd(datacount,p) = max(volTmp(idx_Pulse(cidx)+1:idx_Pulse(cidx)+cdatapts)) - volTmp(idx_Pulse(cidx));  % max voltage, i.e. eta = max(v(2):v(11)) - v(1)
                        currcd(datacount,p) = mean(currTmp(idx_Pulse(cidx)+1:idx_Pulse(cidx)+cdatapts));  % average pulse current
                        % discharge pulses
                        volcd(datacount,p+nPulses) = min(volTmp(idx_Pulse(didx)+1:idx_Pulse(didx)+ddatapts)) - volTmp(idx_Pulse(didx));    % max voltage, i.e. eta = max(v(2):v(11)) - v(1)
                        currcd(datacount,p+nPulses) = mean(currTmp(idx_Pulse(didx)+1:idx_Pulse(didx)+ddatapts));    % average pulse current
                    end

            end % End of HPPC/ SOC set loop
        end % End of if conditon HPPC/SoC start condition
    end % End of cell number loop
end % End of temperature loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Plot current dependency %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
    for k = 1:length(HPPCStart) % loop for each SOC
        subplot(2,2,k);
        
        for j = 1:length(cellNumbers) % loop for each cell
            
            if j==1
                sty='- .'; % cell 1 style
            end
            if j==2
                sty='-- o'; % cell 2 style
            end
            if j==3
                sty=': *'; % cell 3 style
            end
                        
            for i = 1: length(tempValues) % loop for each Temperature 
               
                ri = ((i-1)*length(cellNumbers)+(j-1))*length(HPPCStart)+k; % row index for voltage and current arrays
                
                if i==1
                    plot(currcd(ri,:),volcd(ri,:), sty)
                    legend([num2str(tempValues(i)),'C'])
                    hold on;
                end
                if i==2
                    plot(currcd(ri,:),volcd(ri,:),strcat(sty,'r'))
                    legend([num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
                    hold on;
                end
                if i==3
                    plot(currcd(ri,:),volcd(ri,:),strcat(sty,'g'))
                    legend([num2str(tempValues(i-2)),'C'],[num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
                    hold on;
                end
            end

            xlabel('Current (A)'); ylabel('eta (V)')
            title(['Current depedency SOC: ',num2str(round(socHPPC(k))),'%'])
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Plot Gain %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
    for k = 1:length(HPPCStart) % loop for each SOC
        subplot(2,2,k);
        
        for j = 1:length(cellNumbers) % loop for each cell
            
            if j==1
                sty='- .'; % cell 1 style
            end
            if j==2
                sty='-- o'; % cell 2 style
            end
            if j==3
                sty=': *'; % cell 3 style
            end
                        
            for i = 1: length(tempValues) % loop for each Temperature 
               
                ri = ((i-1)*length(cellNumbers)+(j-1))*length(HPPCStart)+k; % row index for voltage and current arrays
                
                if i==1
                    plot(currcd(ri,:),-volcd(ri,:)./currcd(ri,:), sty)
                    legend([num2str(tempValues(i)),'C'])
                    hold on;
                end
                if i==2
                    plot(currcd(ri,:),-volcd(ri,:)./currcd(ri,:),strcat(sty,'r'))
                    legend([num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
                    hold on;
                end
                if i==3
                    plot(currcd(ri,:),-volcd(ri,:)./currcd(ri,:),strcat(sty,'g'))
                    legend([num2str(tempValues(i-2)),'C'],[num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
                    hold on;
                end
            end

        xlabel('Current (A)'); ylabel('Gain (Ohms)')
        title(['Gain SOC: ',num2str(round(socHPPC(k))),'%'])
        end
    end
    
end

% %old
% %%
% datacount = 1 ; % initialise data counter
% volcd  = NaN([length(tempValues)*4,10]) ; % initialise voltage array for plots
% currcd = NaN([length(tempValues)*4,10]) ; % initialise current array for plots
% for ii = 1: length(tempValues) % loop for each temperature 
% 
%     [dummy,idxTemperature] = ismember(tempValues(ii),tempBP);
%    
%     for jj = 1:length(cellNumbers) % loop for cells NOT USED YET
%      
%         % load data
%         if tempValues(ii)<0
%             tmp = ['m',num2str(abs(tempValues(ii)))];
%             fileLocation = [tmp,'degC'];
%         else
%             fileLocation = [num2str(tempValues(ii)),'degC'];
%         end
%         matFileName = ['Profile_HPPC_Cell',num2str(cellNumbers(jj)),'_',num2str(tempValues(ii)),'degC'];
%         load([fileLocation,'\',matFileName])
%         
%         timeHours = timeVec/3600;   % Convert time vector into hours
%         currVec = -currentVec;      % -ve is charge
%         
%         % Obtain indices of the excitation and relaxation test portion
%         [~,idx_ExcRelax] = find_exc_segments(currVec,1);
%         HPPCStart = hPPCStartMatrix{idxTemperature,cellNumbers(jj)};
%         
%         if ~isempty(HPPCStart)
%             % Obtain cumulative capacity from Bitrode tester to calculate SoC.
%             % Capacity is calulcated from when fully charged to fully discharged
%             for cc = 1:length(startEndRegions)
%                 pp = startEndRegions(cc);
%                 idxTmp = idx_ExcRelax(pp,1):idx_ExcRelax(pp,2);
%                 capTmp = capVec(idxTmp);            % Extract corresponding capacity signal
%                 maxCapDel = max(abs(capTmp));       % Get maximum change in Ah to calculate change in SoC
%                 % Store cummulative charge and discharge capacity for each cell
%                 % at each temperature
%                 if cc == 1
%                     cumCap(cc,jj,ii) = maxCapDel;
%                 else
%                     cumCap(cc,jj,ii) = cumCap(cc-1,jj,ii) + maxCapDel;
%                 end
%             end
%             normCap = Cn; % max(cumCap(:,jj,ii));
%             socEstimate =  (normCap - cumCap(:,jj,ii))./normCap;
%             [~,idx_soc_HPPCStart]= ismember(HPPCStart - 1,startEndRegions);    % Get the indices prior to HPPC cycle start to extract SoC at that point
%             socHPPC(:,jj,ii) = socEstimate(idx_soc_HPPCStart)*100;            % The correponding soc prior to each HPPC cycle for each cell and temperature
%             
%             for kk = 1:length(HPPCStart) % loop for each SOC 
%                 
%                 % Check number of Pulses = 10
%                 if length(idx_ExcRelax)+1 - HPPCStart(kk) < 12
%                     nPulses = floor((length(idx_ExcRelax)+1 - HPPCStart(kk))/2);
%                 else
%                     nPulses = 5;
%                 end
%                
%                 % Select HPPC pulses for a given SoC
%                 ps = HPPCStart(kk);
%                 idxTmp = idx_ExcRelax(ps,1):idx_ExcRelax(ps+2*nPulses-1,2);
%                 timeTmp = timeVec(idxTmp);
%                 currTmp = currVec(idxTmp);
%                 volTmp = volVec(idxTmp);
%                                
% 
%                   
%                     % Obtain indices of the excitation and relaxation test portion
%                     [~,idx_Pulse] = find_exc_segments(currTmp,1);
%                     count = 0 ;                    
%                     for p = 1:nPulses  % loop for each pulse
%                         
%                         cidx = 2*(nPulses+1-p) ; % charge pulse select: 10, 8, 6, 4, 2
%                         didx = 2*p-1 ;           % discharge pulse select: 1, 3, 5, 7, 9
%                         
%                         % check number of sample points on current
%                         % for charge pulses
%                         cdatapts = 1 ; 
%                         while abs(currTmp(idx_Pulse(cidx)+cdatapts+1) - currTmp(idx_Pulse(cidx))) > 1
%                             cdatapts = cdatapts + 1 ;
%                         end
%                         
%                         % for discharge pulses
%                         ddatapts = 1 ; 
%                         while abs(currTmp(idx_Pulse(didx)+ddatapts+1) - currTmp(idx_Pulse(didx))) > 1 
%                             ddatapts = ddatapts + 1 ;
%                         end
% 
%                         % compute max voltage and average current
%                         % charge pulses
%                         volcd(datacount,p) = max(volTmp(idx_Pulse(cidx)+1:idx_Pulse(cidx)+10)) - volTmp(idx_Pulse(cidx));  % max voltage, i.e. eta = max(v(2):v(11)) - v(1)
%                         currcd(datacount,p) = mean(currTmp(idx_Pulse(cidx)+1:idx_Pulse(cidx)+cdatapts));  % average pulse current
%                         % discharge pulses
%                         volcd(datacount,p+nPulses) = min(volTmp(idx_Pulse(didx)+1:idx_Pulse(didx)+10)) - volTmp(idx_Pulse(didx));    % max voltage, i.e. eta = max(v(2):v(11)) - v(1)
%                         currcd(datacount,p+nPulses) = mean(currTmp(idx_Pulse(didx)+1:idx_Pulse(didx)+ddatapts));    % average pulse current
%                     end
%                 
%                 datacount = datacount + 1; % increment data counter
%       
%             end % End of HPPC/ SOC set loop
%         end % End of if conditon HPPC/SoC start condition
%     end % End of cell number loop
% end % End of temperature loop
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% Plot current dependency %%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     figure()
%     for k = 1:length(HPPCStart) % loop for each SOC
%         subplot(2,2,k);
%         
%         for i = 1: length(tempValues) % loop for each Temperature 
%             if i==1
%                 plot(currcd(k,:),volcd(k,:),'- .')
%                 legend([num2str(tempValues(i)),'C'])
%                 hold on;
%             end
%             if i==2
%                 inc = length(HPPCStart) + k ; 
%                 plot(currcd(inc,:),volcd(inc,:),'- .r')
%                 legend([num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
%                 hold on;
%             end
%             if i==3
%                 inc = length(HPPCStart) + inc ;
%                 plot(currcd(inc,:),volcd(inc,:),'- .g')
%                 legend([num2str(tempValues(i-2)),'C'],[num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
%             end
%         end
%         
%         xlabel('Current (A)'); ylabel('eta (V)')
%         title(['Current depedency SOC: ',num2str(round(socHPPC(k))),'%'])
%     end
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% Plot Gain %%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure()
%     for k = 1:length(HPPCStart) % loop for each SOC
%     subplot(2,2,k);
%         
%         for i = 1: length(tempValues) % loop for each Temperature 
%             if i==1
%                 plot(currcd(k,:),volcd(k,:)./currcd(k,:),'- .')
%                 legend([num2str(tempValues(i)),'C'])
%                 hold on;
%             end
%             if i==2
%                 inc = length(HPPCStart) + k ; 
%                 plot(currcd(inc,:),volcd(inc,:)./currcd(inc,:),'- .r')
%                 legend([num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
%                 hold on;
%             end
%             if i==3
%                 inc = length(HPPCStart) + inc ;
%                 plot(currcd(inc,:),volcd(inc,:)./currcd(inc,:),'- .g')
%                 legend([num2str(tempValues(i-2)),'C'],[num2str(tempValues(i-1)),'C'],[num2str(tempValues(i)),'C'])
%             end
%         end
%         
%         xlabel('Current (A)'); ylabel('Gain (Ohms)')
%         title(['Gain SOC: ',num2str(round(socHPPC(k))),'%'])
%     end
%     
% end

