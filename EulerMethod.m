% Numerical methods: Example 1
% Usage w = EulerMethod(@deriFcn,tSpan,w0,N)

function [t,w] = EulerMethod(deriFcn,tSpan,w0,N)

h = (tSpan(2)-tSpan(1))/N;
w(1) = w0;
t(1) = tSpan(1);

for ii = 1:N
    t(ii+1) = tSpan(1)+ ii*h;
    w(ii+1) = w(ii) + h*deriFcn(t(ii),w(ii));
end

