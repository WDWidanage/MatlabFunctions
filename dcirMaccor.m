function [res,options] = dcirMaccor(data,varargin)
% Calculates DCIRs based off MAccor measurements

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments
addRequired(parObj,'data');
addParameter(parObj,'resTimes',[0.4,10]);
addParameter(parObj,'Cn',1);

parse(parObj,data,varargin{:});

data = parObj.Results.data;

resTimes = parObj.Results.resTimes;
timeVec = data.TestTime_s_ - data.TestTime_s_(1);
Cn = parObj.Results.Cn;

% Check if pulse is charge or discharge
[~,idx] = max(abs(data.Amps));
if data.Amps(idx)>0
    options.CurrSgnStr = 'C';
    options.CurrSgn = 1;
else
    options.CurrSgnStr = 'D';
    options.CurrSgn = -1;
end

for rr = 1:length(resTimes)
     if max(timeVec) >= resTimes(rr)
        [~,idx] = min(abs(timeVec - resTimes(rr)));
        res(rr) = abs(data.Volts(idx)-data.Volts(1))/abs(data.Amps(idx));
        
        if (abs(data.Amps(idx))) <  0.95* max(abs(data.Amps))
            options.msg = 'Possible Amps derate';
        else
            options.msg = "Valid resistance";
        end
    else
        res(rr) = nan;
        options.msg = 'Pulse length shorter than request';
    end
    options.current = max(abs(data.Amps));
    options.Crate = max(abs(data.Amps))/Cn; 
    


    
    if round(options.Crate) < 1
        nZ = -floor(log10(options.Crate));
        cInt = round(options.Crate*10^nZ);
        if options.current == 0
            options.Str = 'Pulse0C';
        else
            options.Str = ['PulseC',repmat('0',1,nZ-1),num2str(cInt)];
        end
    else        
        options.Str = ['Pulse',num2str(round(options.Crate)),'C'];
    end
end

end