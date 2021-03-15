function [GJW, PO, M, NoiseSD ExpYL,ExpXL]=EFRF(ysig,xsig,L,Harmonics,wind,Bd)
%Calculates the estimated transfer function as the ratio of cross power
%spectrum and auto power spectrum, for different amounts of overlap. WOSA
%method

%Calculates the estimated transfer function as the ratio of cross power
%spectrum and auto power spectrum, for different amounts of overlap.

%Bd: Block displacment. If null, all possible block displacments are
%performed.
%wind: type of window applied
try dummy=Bd; catch Bd=[]; end;

Nsig=length(xsig)/L; %Number of blocks in the signal without any overlap
NOS=find(mod(L,[1:L/2])==0);%nonoverlap samples
NOS=fliplr(NOS);
NO=round((Nsig-1)*L./NOS)+1;%Number of blocks with ovelapping, NO~=Nsig
nb=[Nsig:NO(1)-1]; %Number of possible blocks between 0% and 50% overlap, 
%However last block does not coincide with end of record
%Calculating the number of non-overlapping samples for these blocks 
ns=round((Nsig-1)*L./(nb-1));
NOB=[nb,NO]; %Augment the number of blocks
if Bd
   NOS=Bd; 
   NOB=round((Nsig-1)*L./NOS)+1;%Number of blocks with ovelapping, NO~=Nsig
else
   NOS=[ns,NOS]; %Augment the number of non overlapping samples
end

PO=100*(L-NOS)/L;

for i=1:length(NOS)%Loop for each percentage overlap
     nos=NOS(i);nob=NOB(i);
     SYX=zeros(1,length(Harmonics));SXX=zeros(1,length(Harmonics));t2=0;k=0;%Intitialise for loop
     EYL=zeros(1,length(Harmonics));SYL=zeros(1,length(Harmonics));
     EXL=zeros(1,length(Harmonics));
     %now treat problem as normal FRF problem and evaluate over NO blocks
     while (t2+L-1)<=length(xsig)
         k=k+1;
         t1=round((k-1)*nos)+1;t2=t1+L-1;
         yL=ysig(t1:t2).*wind;xL=xsig(t1:t2).*wind;% Take a block of length L and window it
         YL=fft(yL);XL=fft(xL);%FFT with window
         SYX=SYX+YL(Harmonics+1).*conj(XL(Harmonics+1));%Cross power spectrum
         SXX=SXX+XL(Harmonics+1).*conj(XL(Harmonics+1));%Auto power spectrum
         %Calculating complex variance from each YL block at the specified
         %harmonics
         EYL=EYL+(YL(Harmonics+1))/nob;SYL=SYL+(YL(Harmonics+1).*conj(YL(Harmonics+1)))/nob;
         EXL=EXL+(XL(Harmonics+1))/nob;
         t2=t1+nos;%update t2 to the next position of t1
     end
    M(i)=k;
    GJW(i,:)=SYX./SXX; %Estimated FRF
    NoiseSD(i,:)=sqrt(SYL-(EYL.*conj(EYL)));%Estimated noise standard deviation
    ExpYL(i,:)=EYL;
    ExpXL(i,:)=EXL;
end

%ALL COMENTED FROM PREVIOUS CODING
% for i=1:length(NOS)
%     nos=NOS(i);
%     NO=round((Nsig-1)*L/nos)+1;%Number of blocks with ovelapping, NO~=Nsig
%     SYX=zeros(1,L);SXX=zeros(1,L);
%     EYL=zeros(1,length(Harmonics));SYL=zeros(1,length(Harmonics));
%     EXL=zeros(1,length(Harmonics));
%     %now treat problem as normal FRF problem and evaluate over NO blocks
%     for k=1:NO
%         t1=round((k-1)*nos)+1;t2=t1+L-1;
%         yL=ysig(t1:t2);xL=xsig(t1:t2);%Take a block of L sample
%         YL=fft(yL);XL=fft(xL);%FFT without any window
%         SYX=SYX+YL.*conj(XL);%Cross power spectrum
%         SXX=SXX+XL.*conj(XL);%Auto power spectrum
%         %Calculating complex variance from each YL block at the specified
%         %harmonics
%         EYL=EYL+(YL(Harmonics+1))/NO;SYL=SYL+(YL(Harmonics+1).*conj(YL(Harmonics+1)))/NO;
%         EXL=EXL+(XL(Harmonics+1))/NO;
%     end
%     GJW(i,:)=SYX./SXX; %Estimated FRF
%     NoiseSD(i,:)=sqrt(SYL-(EYL.*conj(EYL)));%Estimated noise standard deviation
%     ExpYL(i,:)=EYL;
%     ExpXL(i,:)=EXL;
% end