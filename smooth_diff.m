function dy = smooth_diff(x,y,n,g)
% SMOOTH_DIFF  Derivative of noisy data.
% This function takes a set of noisy data (x,y) and calculates dy/dx. The
% technique used is fitting a degree p polynomial tothe n neighbours to 
% each side of each data point and taking the polynomial to calculate the
% derivative. It returns dy which is the derivative
%   dy = SMOOTH_DIFF(x,y,n,p) returns the derivative of the data set (x,y)
%   fitting a polynomial of degree p to the n neighbours on each side.

dy = zeros(size(y));

p = polyfit(x(1:2*n+1),y(1:2*n+1),g);
dp = polyder(p);
dy(1:n+1) = polyval(dp,x(1:n+1));

for i = (n+2):(length(y)-n-2)
    p = polyfit(x(i-n:i+n),y(i-n:i+n),g);
    dp = polyder(p);
    dy(i) = polyval(dp,x(i));
end

p = polyfit(x(end-2*n:end),y(end-2*n:end),g);
dp = polyder(p);
dy(end-n-1:end) = polyval(dp,x(end-n-1:end));

    
end