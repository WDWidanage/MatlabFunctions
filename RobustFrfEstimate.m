function varargout=RobustFrfEstimate(U,Y)
%
%Robust method for estimating a frequency response. Input data are periodic
%measured over several periods and repeated over several realisations
%(experiments). 
%
%Input arguments
%
%   U = FFT of input signals, size F x P x R. Rows are the frequencies at which
%   the frf is estimated. Columns are the periods. Finial dimension for
%   realisations.
%   Y = FFT of output signal, size F x P x R similar to U
%
%Output arguments
%   Gfrf=RobustFrfEstimate(U,Y)
%   [Gfrf,stdG]=RobustFrfEstimate(U,Y)
%   [Gfrf,stdG,stdGn]=RobustFrfEstimate(U,Y)
%   [Gfrf,stdG,stdGn,stdGs]=RobustFrfEstimate(U,Y)
%
%   Gfrf  = Estimated frf, best linear approximation, size F x 1
%   stdG  = Standard deviation of the frf, total variance, size F x 1
%   stdGn = noise standard deviation on the frf, size F x 1
%   stdGs = stocahstic nonlinear distortion variance, size F x1%
%
%W.D.Widanage (01/05/2010) (@home)

[F,P,R]=size(U);

G_EveryPeriod=Y./U;                         %frf for each period and each realisation
Gfrf=mean(mean(G_EveryPeriod,2),3);         %best linear approximate

varG=var(mean(G_EveryPeriod,2),0,3)/R;      %variance over averaged periods scaled by R to give
                                            %variance of Gfrf which is a statistic 
                                            %averaged over realisations
stdG=sqrt(varG);
                                        
varGn=mean(var(G_EveryPeriod,0,2),3)/(P*R); %variance over periods and the mean over realisations,
                                            %scaled by P and R to make it
                                            %relative to the statistic Gfrf
stdGn=sqrt(varGn);

varGs=R*(varG-varGn);                       % variance of stocastic nonlinear contributions
stdGs=sqrt(varGs);

switch nargout
    case 1
        varargout={Gfrf};
    case 2
        varargout={Gfrf,stdG};
    case 3
        varargout={Gfrf,stdG,stdGn};
    case 4
        varargout={Gfrf,stdG,stdGn,stdGs};
end
end