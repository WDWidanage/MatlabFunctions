
function u=SkewMultisine(U,order,desired,N,RMS)

% Designing a multisine with a desired amplitude distribution
%
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
%
% W.D.Widanage 09/07/10 (Bible Black)
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


if A==0
    if U(1)==0
        A=sqrt(N^2*RMS^2/(2*(lU-1))); %Amplitudes
    else
        A=sqrt(N^2*RMS^2/(2*lU)); %Amplitudes
    end
end

udes = desired;

%Amplitude vector
Amp(U+1)=A;
Amp(N-U+1)=A;
if U(1)==0
    Amp(end)=[];
end

if mod(N,2)==0
    Nq=N/2;
else
    Nq=(N+1)/2-1;
end

%Regressor matrix
K=zeros(N,N);
if mod(N,2)==0;
    for nn=0:N-1
        K(nn+1,:)=[Amp(1),Amp(2:Nq).*exp(j*2*pi/N*nn*[1:Nq-1])+Amp(2:Nq).*exp(j*2*pi/N*nn*[N-1:-1:Nq+1]),0,j*Amp(2:Nq).*exp(j*2*pi/N*nn*[1:Nq-1])-j*Amp(2:Nq).*exp(j*2*pi/N*nn*[N-1:-1:Nq+1])];
    end
    zeroCols = find(all(K==0));
    K(:,zeroCols)=[];
else
    for nn=0:N-1
        K(nn+1,:)=[Amp(1),Amp(2:Nq+1).*exp(j*2*pi/N*nn*[1:Nq])+Amp(2:Nq+1).*exp(j*2*pi/N*nn*[N-1:-1:N-Nq]),j*Amp(2:Nq+1).*exp(j*2*pi/N*nn*[1:Nq])-j*Amp(2:Nq+1).*exp(j*2*pi/N*nn*[N-1:-1:N-Nq])];
    end
    zeroCols = find(all(K==0));
    K(:,zeroCols)=[];
end


%Least squares approx
expTheta=real(K\udes);

ComPar=[expTheta(1:max(U))+j*expTheta(max(U)+1:end)];

%Normalise magnitude
NormPhase=ComPar;%./abs(ComPar);

%Generate new time signal
%Frequency components
Au=NormPhase;

Vtemp=zeros(N,1);

%Frequency signals
Vtemp(U+1)=Au;
Vtemp(N-U+1)=conj(Au);
if U(1)==0
    Vtemp(end)=[];
end

%Initial time signal
u=ifft(Vtemp)*N;
end