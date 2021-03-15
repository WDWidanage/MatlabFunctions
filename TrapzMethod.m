% Numerical methods: Example 1
% Usage w = TrapzMethod(@deriFcn,tSpan,w0,N)

function [t,w] = TrapzMethod(deriFcn,tSpan,w0,N)

h = (tSpan(2)-tSpan(1))/N;
w(1) = w0;
t(1) = tSpan(1);

for ii = 1:N
    s1 = deriFcn(t(ii),w(ii));
    t(ii+1) = tSpan(1)+ ii*h;
   
    wt = w(ii) + h*s1;              % One-step Euler prediction
    s2 = deriFcn(t(ii+1),wt);
    w(ii+1) = w(ii) + h/2*(s1+s2);
end

