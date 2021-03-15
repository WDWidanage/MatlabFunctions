function [NpClk,harmVec,timeVec,freqVec] = GenMultiLevN(Tp,varargin)

% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 12/08/2019 (Shimmer)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only


p = inputParser; % Create an input parse object to handle positional and property-value arguments


% Create variable names and assign default values after checking the value
addRequired(p,'Tp', @isnumeric);            % Signal period in hours

% Optional parameters
addParameter(p,'clkTs',1)                  % Default clock interval in seconds: 1s
addParameter(p,'fMax',2E-3)                % Default maximum frequency in mHz: 2mHz
addParameter(p,'suppHarmMul',0)            % Default multiples of suppressed haromonics: None
addParameter(p,'Ts',1)                     % Defualt sampling time in seconds: 1s

% Re-parse parObj
parse(p,Tp, varargin{:})
Tp = p.Results.Tp;
clkTs = p.Results.clkTs;
Ts = p.Results.Ts;
suppHarmMul = p.Results.suppHarmMul;
fMax = p.Results.fMax;


NTmp = floor(Tp*3600/clkTs);        % Number of clock pulses
hMinClk = 1;

if suppHarmMul == 0
    NpClk = NTmp;
else
    prodHarmMul = prod(suppHarmMul);
    rem = mod(NTmp,prodHarmMul);
    NpClk = NTmp - rem;
end

f0 = 1/(NpClk*clkTs);                      % Fundemental frequency
hMaxClk = floor(fMax/f0);


if suppHarmMul == 0
    harmVec = [hMinClk:hMaxClk]';
else
    harmVecTmp = [hMinClk:hMaxClk]';
    for ii = 1:length(suppHarmMul)
        hM = suppHarmMul(ii);
        harmVecTmp = setdiff(harmVecTmp,[hM:hM:hMaxClk]);
        harmVec = harmVecTmp';
    end
end
if NpClk < 2*hMaxClk 
    error('Maximum specified frequency must be less than %2.2E Hz',NpClk*f0/2); % Np/(2*Tp*3600)
end
if mod(clkTs,Ts) ~= 0 
    tStep = clkTs - mod(clkTs,Ts);
    error('Sampling time (Ts) is not a multiple of timeStep.\nNearest multiple for timeStep: %3ds',tStep); 
end

timeVec = [0:NpClk-1]'*clkTs;
freqVec = harmVec/(NpClk*clkTs);

end

