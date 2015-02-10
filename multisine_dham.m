
function u=multisine_dham(U,order,N,RMS)

%Output 
%u: Samples of the time signal.
%
%Mandatory arguments
%U: Vector of the input harmonics eg U=[1:2:50]
%order: a scaler value that defines the Nyquist frequency point as
%       order x highest harmonic number in U. Thus creating a margin
%       allowing nolinearties upto that order to be detected at the output
%       signal, eg order=1. As this defines the Nyquist frequency the
%       minimum signal length Nmin is also determined.
%
%Optional arguments
%N: To define a period different from that of Nmin The value has to be at
%   least twice the highest harmonic number specifed in U. N can also be
%   set to zero treating as if no period is specified and allowing the 
%   period length to be determeined by the order and a desired rms value 
%   may be set. 
%RMS: The desired root mean square value of the signal. If the signal has a
%     bias, this is not considered and it is the rms about the bias value.
%     This is also the standard deviation of the signal.

try 
    N;
    SP=0;
    if N==0; %Treat as if signal period is not passed as an input, SP=1 to specify period later on 
        SP=1; 
    end
catch
    SP=1; %If signal period is not passed as an input, SP=1 to specify period later on
end;
try
    RMS; A=0; catch A=1; %If root mean square is not passed amplitudes=1
end


%Change to column vectors
U=U(:);
lU=length(U);

%Highest harmonic
Hmax=max(U);

if SP==1
    %Nyquist frequency point inorder to detect the higest harmonic of the
    %output.
    Nf=order*Hmax;
    %Therefore required smallest signal period N
    N=2*Nf+1;
end

Utemp=zeros(N,1);

%Phases
if U(1)==0
    Thetatemp1=0;
    Thetatemp2=2*pi*rand(lU-1,1);
    ThetaU=[Thetatemp1;Thetatemp2];
else
    ThetaU=2*pi*rand(lU,1);
end

if A==0
    if U(1)==0
        A=sqrt(N^2*RMS^2/(2*(lU-1))); %Amplitudes
    else
        A=sqrt(N^2*RMS^2/(2*lU)); %Amplitudes
    end
end
%Frequency components
Au=A*exp(j*ThetaU);


%Frequency signals
Utemp(U+1)=Au;
Utemp(N-U+1)=conj(Au);
if U(1)==0
    Utemp(end)=[];
end

%Time signals
u=ifft(Utemp);

end
