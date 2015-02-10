function [Model]=San_algorithm(G,w,nb,na,A,domain,fs,T)

%Check if domain, sampling frequency and number of iterations are provided
try 
    domain;
catch
    domain='z';
end
try 
    fs;
catch
    fs =1;
end
try 
    T;
catch
    T =20;
end

j=sqrt(-1);
if domain=='z'
    v=exp(-j*w);    
else
     v=j*w*fs;
end


F=length(w);
%Construct regressor matrix
for kk=1:F
    Reg(kk,:)=[-G(kk)*(v(kk).^(na:-1:1)) v(kk).^(nb:-1:0)];
end


for ll=1:T
    A_=diag(1./A);
    RegA=A_*Reg;GA=A_*G;

    Regt=[real(RegA);imag(RegA)];
    Gt=[real(GA);imag(GA)];

    %Numerically stabe least squares
    Scale=sqrt(sum(Regt.^2)); %l2 norm of each column
    idxZeros=find(Scale<1E-14);
    Scale(idxZeros)=1;
    Regt=Regt./repmat(Scale,2*F,1);
    %SVD of Reg
    [Un, Sn, Vn] = svd(Regt,0);
    ss = diag(Sn);
    idxZeros = find(ss < 1E-14);
    ss(idxZeros) = inf;
    ss = diag(1./ss);
    theta=Vn*ss*Un'*Gt;
    theta=theta./Scale';

    A=polyval([theta(1:na);1],v);
end

Model.A=theta(1:na);
Model.B=theta(na+1:end);