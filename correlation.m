function [R,tau]=correlation(x,y,R1,R2,plotg)

%Computes the auto or cross corelation function for tau between R1 and R2
%x,y : column vectors of input and output, of same length, real or complex
%R1,R2: lower and upper limit for tau
%plotg: if plotg is not a null variable, plots correlation vs tau 
%R: Correlation matrix of length R2-R1+1
%tau: vector of tau

try dummy=plotg; catch plotg=[]; end;

L=length(x);
for r=0:(R2-R1)
    tlen=L-abs(r+R1);
    if r<=abs(R1)
        R(r+1)=sum(y(1:tlen).*x(L-tlen+1:L))/(tlen+1);
%         R(r+1)=sum(conj(y(1:tlen)).*x(L-tlen+1:L))/L;
    else 
        R(r+1)=sum(y(L-tlen+1:L).*x(1:tlen))/(tlen+1);
%         R(r+1)=sum(conj(y(L-tlen+1:L)).*x(1:tlen))/L;
    end
end
tau=[R1:R2];

if plotg
    figure(plotg);
    if x==conj(x), plot(tau,R);
    xlabel('Time shift'), ylabel('Correlation function');grid on; 
    else plot(R,'+'); 
    xlabel('Real'), ylabel('Imaginary');grid on;
    end;
    
end