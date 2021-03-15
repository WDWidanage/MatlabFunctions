function [Model] = levi_algorithm(G,w,nb,na, varargin)
%Check if domain and sampling frequency are provided

p = inputParser;

% Create variable names and assign default values after checking the value
addRequired(p,'G', @isnumeric);
addRequired(p,'w', @isnumeric);
addRequired(p,'nb', @isnumeric);
addRequired(p,'na', @isnumeric);


% Optional parameteres
addParameter(p,'domain','s');
addParameter(p,'fs',1);


% Re-parse parObj
parse(p,G,w,nb,na,varargin{:})


j=sqrt(-1);
if domain=='s' 
    v = j*w*fs;
else
    v = exp(-j*w);
end


F=length(w);

%Construct regressor matrix
for kk=1:F
    Reg(kk,:) = [-G(kk)*(v(kk).^(na:-1:1)) v(kk).^(nb:-1:0)];
end

Regt = [real(Reg);imag(Reg)];
Gt = [real(G);imag(G)];

theta = Lls(Regt,Gt);

%Numerically stabe least squares
% Scale=sqrt(sum(Regt.^2)); %l2 norm of each column
% idxZeros = Scale<1E-14;
% Scale(idxZeros)=1;
% Regt=Regt./repmat(Scale,2*F,1);
% 
% %SVD of Reg
% [Un, Sn, Vn] = svd(Regt,0);
% ss = diag(Sn);
% idxZeros = ss < 1E-14;
% ss(idxZeros) = inf;
% ss = diag(1./ss);
% theta=Vn*ss*Un'*Gt;
% theta=theta./Scale';


Model.A = theta(1:na);
Model.B = theta(na+1:end);