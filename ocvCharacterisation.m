function ocvResults = ocvCharacterisation(currVec, volVec, timeVec,varargin)
% Get OCV vs capacity (Ah) and SoC
%
% Copyright (C)  W.D. Widanage   -  WMG, University of Warwick, U.K. 27/10/19 (Clair de Lune)
% All Rights Reserved



% Create an input parse object to handle positional and property-value arguments
parObj = inputParser;

addRequired(parObj,'currVec')
addRequired(parObj,'volVec')
addRequired(parObj,'timeVec')

addParameter(parObj,'capVec',[]);
addParameter(parObj,'plotSeg',0);
addParameter(parObj,'pulseStrt',3);


% Re-parse parObj - 09-July-2015
parse(parObj,currVec, volVec, timeVec,varargin{:});

currVec = parObj.Results.currVec;
volVec = parObj.Results.volVec;
timeVec = parObj.Results.timeVec;
capVec = parObj.Results.capVec;
plotSeg = parObj.Results.plotSeg;
pulseStrt = parObj.Results.pulseStrt;
timeHours = timeVec/3600;
[~,idx_ExcRelax] = find_exc_segments(currVec,'tol',0.05,'plotSeg',plotSeg,'timeVec',timeHours);


lExcRelax = size(idx_ExcRelax);
numPulses = lExcRelax - pulseStrt + 1;

Qn = 0; % Initialise capacity variable
Qc = 0; % Initialise charge capacity variable
Qd = 0; % Initialise discharge capacity variable
OCV = volVec(idx_ExcRelax(pulseStrt,1)); % Initialise OCV
cntrc = 0;
cntrd = 0;

for pp = 1: numPulses
    idxStrt = idx_ExcRelax(pulseStrt+pp-1,1);
    idxEnd = idx_ExcRelax(pulseStrt+pp-1,2);
    idxTmp = idxStrt:idxEnd;
    currTmp = currVec(idxTmp);
    timeTmp = timeHours(idxTmp);
    OCV(pp+1,1) = volVec(idxEnd);
    
    % Compute capacity in Ah for applied pulse
    if isempty(capVec)
        [~,capTmp] = CoulombCounting(currTmp,timeTmp,100,100,0);
    else
        capTmp = capVec(idxTmp);
    end
    
    % Find sgn of current +ve for charge -ve for discharge
    currSgn = currSgnFcn(currTmp);
    
    % Get maximum change in Ah to calculate change in SoC
    maxCapDel = max(abs(capTmp));
    
    % Store cummulative charge and discharge capacity
    
    if currSgn > 0                          % A charge pulse
        cntrc = cntrc + 1;
        if cntrd > 1                        % If any discharge pulses have occured
            Qn(pp+1,1) = Qn(pp) - maxCapDel;
            if pp-cntrd == 1
                OCVc(1) = volVec(idxStrt);
                Qc(pp-cntrd) = Qd(end);
            end
            Qc(pp+1-cntrd ,1) = Qc(pp-cntrd) - maxCapDel;
            OCVc(pp+1-cntrd,1) = volVec(idxEnd);
        else
            Qn(pp+1,1) = Qn(pp) + maxCapDel;
            if pp == 1
                OCVc(1) = volVec(idxStrt);
            end
            Qc(pp+1,1) = Qc(pp) + maxCapDel;
            OCVc(pp+1,1) = volVec(idxEnd);
        end
    else                                    % A discharge pulse
        cntrd = cntrd + 1;
        if cntrc > 1                        % If any charge pulses have occured
            Qn(pp+1,1) = Qn(pp) - maxCapDel;
            if pp-cntrc == 1
                OCVd(1) = volVec(idxStrt);
                Qd(pp-cntrc) = Qc(end);
            end
            Qd(pp+1-cntrc,1) = Qd(pp-cntrc) - maxCapDel;
            OCVd(pp+1-cntrc,1) = volVec(idxEnd);
        else
            Qn(pp+1,1) = Qn(pp) + maxCapDel;
            if pp == 1
                OCVd(1) = volVec(idxStrt);
            end
            Qd(pp+1,1) = Qd(pp) + maxCapDel;
            OCVd(pp+1,1) = volVec(idxEnd);
        end
    end
    
end
% Perform interpolations
refSoC = [0:100]';
Cn = max(Qd);
socC = (Cn-Qc)/Cn*100;
socD = (Cn-Qd)/Cn*100;
refOCVd = interp1(socD,OCVd,refSoC);
refOCVc = interp1(socC,OCVc,refSoC);

% Mean OCV
meanRefOCV = mean([refOCVd,refOCVc],2);

ocvResults.Qn = Qn;
ocvResults.Qc = Qc;
ocvResults.Qd = Qd;
ocvResults.OCV = OCV;
ocvResults.OCVc = OCVc;
ocvResults.OCVd = OCVd;
ocvResults.refOCVd = refOCVd;
ocvResults.refOCVc = refOCVc;
ocvResults.meanRefOCV = meanRefOCV;
ocvResults.refSoC = refSoC;
ocvResults.hystersis = [refOCVc-refOCVd];



end

function currSgn = currSgnFcn(currVec)
[~,idxMax] = max(abs(currVec));
currSgn = sign(currVec(idxMax));
end

function [idx_Exc,idx_ExcRelax] = find_exc_segments(current,varargin)
% This function is for cell pulse testing. Based on a current load profile the
% function returns the index pairs of the excitation segments and
% excitation with relaxation segements. The function provides all the
% excitation and relaxation segments witin the first and last rest
% intervals
%
% Mandotory input arguments:
%   current: Current vector (A), vector size N x 1
%
% Optional input arguments. Create a structure variable with the following fields:
%   tol: tolorance level to define zero level, default 0.1A, size 1 x 1
%   plotSeg: set to 1 to plot ansd display signal segments, default 0, size 1 x 1
%   timeVec: time column of measured current. This is used for plotting, vector size N x 1
%
% Output arguments:
%   idx_Exc: Index pairs for excitation segment, matrix p x 2
%   idx_ExcRelax: Index pairs for excitation segment with relaxation, matrix p x 2
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 24/09/2013 (Welcome to the jungle!!)
% All Rights Reserved

pObj = inputParser; % Create an input parse object to handle positional and property-value arguments


% Create variable names and assign default values after checking the value
addRequired(pObj,'current', @isnumeric);

% Optional parameters
addParameter(pObj,'tol',0.1,@isnumeric);
addParameter(pObj,'plotSeg',0,@isnumeric);
addParameter(pObj,'timeVec',[],@isnumeric);


% Re-parse parObj
parse(pObj,current,varargin{:})

current = pObj.Results.current;
tol = pObj.Results.tol;
plotSeg = pObj.Results.plotSeg;
timeVec = pObj.Results.timeVec;

logical_segments = -tol < current & current < tol; % logic 1 => zero segments, logic 0 => excitation segments

% Initialise
if logical_segments(1) == true
    s(1) = 1;       % initialise start index to 1 and set start vector index to 1
    jj = 1;
else
    jj = 0;        % else set start vector index to 0 since the value will arrive later in loop
end

kk = 0;            % set end vector index to 0 since the value will arrive later in loop
for ii = 2:length(current)
    logical_prev = logical_segments(ii-1);
    if logical_segments(ii) ~= logical_prev    % Detect a change in event
        if logical_segments(ii) == true        % True indicates start of zero segment
            jj = jj+1;
            s(jj) = ii;                        % Save index of start segment
        else                                   % False indicates end of zero segment
            kk = kk+1;
            e(kk) = ii-1;                      % Save index of end segment
        end
    end
end

% If last segment is a zero segmment there wont be a change in event in the for loop
% and manually set the last element of the end vector to the signal length
if logical_segments(end) == true
    e(kk+1) = length(current);
end

idx_Exc = [e(1:end-1)',s(2:end)'-1];
idx_ExcRelax = [e(1:end-1)',e(2:end)'];

[lEx, ~] = size(idx_Exc);
[lExR, ~] = size(idx_ExcRelax);

if  isempty(timeVec)
    timeVec = [1:length(current)]';
end

if plotSeg
    figure
    plot(timeVec',current,'. -', timeVec(idx_Exc(:,1)),zeros(lEx,1),'g o','MarkerFaceColor','g');
    hold on;
    plot(timeVec(idx_Exc(:,2)),current(idx_Exc(:,2)),'r o','MarkerFaceColor','r');
    plot(timeVec(idx_ExcRelax(:,2)),zeros(lExR,1),'k o');
    xlabel('Time'); ylabel('Curent')
    hold off;
end

end

function [soc, ampSec] = CoulombCounting(I,time,Cn_c,Cn_d,SoC0)
% Coulomb counting to estimate state-of-charge
% Perform current integration and normalise with resepect to battery
% capacity. Cn_c and Cn_d are the measured battery capacity when charging
% and discharging respectively at a given temperature.
%
% Integration is performed using trapezoidal method with saturation limits
% of 0 and 1
%
% Input arguments:
%   I: Current vector (A), +ve is assumed discahrging, size N x 1
%   time: time vector (s), size N x 1
%   Cn_c: Cell capacity when charging (Ah), size 1 x 1
%   Cn_d: Cell capacity when discharging (Ah), size 1 x 1
%   SoC0: Initial state-of-charge (%), size 1 x 1
%
% Output arguments:
% soc: Remaining state-of-charge (%), size N x 1
% ampSec: Ampere seconds (As), size N x 1
%
% W.D. Widanage 09/02/2013 (Boogie woogie!)

N = length(I);

% Initialise
soc = zeros(N,1);
ampSec = zeros(N,1);
soc(1) = SoC0;
ampSec(1) = 0;
inc = 0;

for ii = 1:N-1
    
    delT = time(ii+1)-time(ii);
    if I(ii)>=0
        inc = ((I(ii)+I(ii+1))*delT/2)/(Cn_d*3600);  % trapezoid increment approximation method and normalise with capacity (As)
    elseif I(ii)<0
        inc = ((I(ii)+I(ii+1))*delT/2)/(Cn_c*3600);  % trapezoid increment approximation method and normalise with capacity (As)
    end
    
    % Perform integration with saturation limits of 1 and 0
    if  soc(ii)-inc >1
        soc(ii+1) = 1;
    elseif soc(ii)-inc <0
        soc(ii+1) = 0;
    else
        soc(ii+1) = soc(ii)-inc;    % Accumulate remaining capacity
    end
    ampSec(ii+1) = ampSec(ii) + (I(ii)+I(ii+1))*delT/2;
end
end


