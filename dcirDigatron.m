function [res,options] = dcirDigatron(data,varargin)
% Calculates DCIRs based off Digatron measurements

parObj = inputParser; % Create an input parse object to handle positional and property-value arguments
addRequired(parObj,'data');
addParameter(parObj,'resTimes',[0.4,10]);
addParameter(parObj,'Cn',1);

parse(parObj,data,varargin{:});

data = parObj.Results.data;

resTimes = parObj.Results.resTimes;
timeVec = data.ProgTime - data.ProgTime(1);
Cn = parObj.Results.Cn;

% Check if pulse is charge or discharge
[~,idx] = max(abs(data.Current));
if data.Current(idx)>0
    options.CurrSgnStr = 'C';
    options.CurrSgn = 1;
else
    options.CurrSgnStr = 'D';
    options.CurrSgn = -1;
end

for rr = 1:length(resTimes)
%     if max(timeVec) >= resTimes(rr)
        [~,idx] = min(abs(timeVec - resTimes(rr)));
        res(rr) = abs(data.Voltage(idx)-data.Voltage(1))/abs(data.Current(idx));
        
        if (abs(data.Current(idx))) <  0.95* max(abs(data.Current))
            options.msg = 'Possible current derate';
        else
            options.msg = 'Valid resistance';
        end
%     else
%         res(rr) = nan;
%         options.msg = 'Pulse length shorter than request';
%     end
    options.current = max(abs(data.Current));
    options.Crate = max(abs(data.Current))/Cn; 
    


    
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