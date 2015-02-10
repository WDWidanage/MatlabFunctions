function Model = GN_algorithm(G,w,B,A,domain,fs)

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

j=sqrt(-1);
if domain=='z'
    v=exp(-j*w);    
else
     v=j*w*fs;
end

F=length(w);
na=length(A);
nb=length(B)-1;
theta=[A;B];

iter=0;
V=1;
while V > 10^-6 & iter<10
    iter=iter+1;
    Bk=polyval(theta(na+1:end),v);
    Ak=polyval([theta(1:na);1],v);
    Gm=Bk./Ak; Gm=Gm(:);
    error=G-Gm;
    V(iter)=error'*error;
    %Construct Jacobian
    for kk=1:F
        Bkk=Bk(kk);
        Akk=Ak(kk);
        dGda=(-Bkk/Akk^2)*(v(kk).^(na:-1:1));
        dGdb=(1/Akk)*(v(kk).^(nb:-1:0));
        J(kk,:)=[dGda,dGdb];
    end
    Regt=[real(J);imag(J)];
    et=[real(error);imag(error)];
    
    %Numerically stable least squares
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
    delta=Vn*ss*Un'*et;
    delta=delta./Scale';
    theta=theta+delta;
end

Model.A=theta(1:na);
Model.B=theta(na+1:end);
