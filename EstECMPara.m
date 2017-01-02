function [thetaOpt, resNorm] = EstECMPara(t,u,v,order,options,plotFit)
% Estimate the parameters of an ECM
%
% W.D. Widanage 07/01/2016 (Fade to Black)


% Ensure inputs are column vectors
t = t(:);
u = -u(:); % Assumes negative current as charging
v = v(:);

% Fit equivalent circuit on time data
OCV = v(1);
dtTmp = diff(t);
dt = [dtTmp(1);dtTmp];
fh = @(theta,u)EquivalentCircuitModel(theta,u,OCV,dt,order);    % Model function handle


% Set initial values for ECM parameters
[~, idx] = max(abs(u));
if sign(u(idx)) < 0                                 % For a charge pulse
    Ro = (max(v) - OCV)/abs(u(idx));
    OCVacc0 = (v(end) - OCV)/(abs(u(idx))*10); % OCV', Steady state voltage = OCVacc*Area under applied current + OCV
else                                                % For a discharge pulse
    Ro = (OCV - min(v))/abs(u(idx));
    OCVacc0 = (OCV - v(end))/(abs(u(idx))*10); % OCV', Steady state voltage = OCVacc*Area under applied current + OCV
end
                        
for mm = 1: order
    Rp(mm) = 1e-3;
    tau(mm) = 10^(mm-1);
end

theta0 = [Ro OCVacc0 Rp tau];
thetaLb = zeros(2*order+2,1);    % Define parameter lower bound to be zero
thetaUb = [];    % Parameter upper bound is undefined
[thetaOpt,resNorm] = lsqcurvefit(fh,theta0,u,v,thetaLb,thetaUb,options);

% Simulate ECM with estimation data set and estimated optimum parameteres
yECM = EquivalentCircuitModel(thetaOpt,u,OCV,dt,order);
errorECM = v - yECM;
rmsErrECM = rms(errorECM);
pkErrECM = max(abs(errorECM));

if plotFit == 1
    figure();
    subplot(2,1,1)
    plot(t,u,'. -');
    ylabel('Current (A)');
    subplot(2,1,2)
    plot(t,v,'- .',t,yECM,'-')
    xlabel('Time (s)'); ylabel('Voltage (V)'); legend('Measured','ECM sim')
    title(['RMSE ECM: ',num2str(rmsErrECM),' Pk Err ECM: ', num2str(pkErrECM)]);
end



function [Vl, J]= EquivalentCircuitModel(theta, Il, OCV, del, order)
%
%
% Parameters are arranged as:
%      theta = Ro, OCV', Rp1,..,Rpn, tau1,...,taun
%
% W. D. Widanage 19/10/2013 (Regain)
%

% Change parameter vector to a column vector
theta = theta(:);

% Extract parameters
Ro = theta(1);
OCVacc = theta(2);
RpAll = theta(3:2+order);
tauAll = theta(3+order:2+2*order);

dataPts = length(Il);       % Number of data points

% Initialise
Ip = zeros(order,1);
Vl = zeros(dataPts,1);
Integral = 0;

% Initialise for Jacobian
J = zeros(dataPts,length(theta));
dIpdtau = zeros(order,1);

Vl(1) = OCV - Ro*Il(1) - OCVacc*Integral - RpAll'*Ip;
J(1,:) = [-Il(1), -Integral, -Ip', -dIpdtau'];

for ii = 2:dataPts
    Ts = del(ii-1);                                    % Sampling interval
    expTau = exp(-Ts./tauAll);
    for jj = 1:order
        Ip(jj,1) = expTau(jj)*(Ip(jj,1) - Il(ii)) + Il(ii);
    end
    Integral = Integral + (Il(ii) + Il(ii-1))*Ts/2;
    
    Vl(ii) = OCV - Ro*Il(ii) - OCVacc*Integral - RpAll'*Ip;
    
    % Create Jacobian matrix
    for jj = 1:order
        dIpdtau(jj,1) = expTau(jj)*(dIpdtau(jj,1) + Ip(jj,1)*Ts/tauAll(jj)^2 - Il(ii)*Ts/tauAll(jj)^2);
    end
    RdIpdtau = RpAll.*dIpdtau;
    J(ii,:) = [-Il(ii), -Integral, -Ip', -RdIpdtau'];
end