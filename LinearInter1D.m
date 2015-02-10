function Vx = LinearInter1D(xBP,V,x)

% Perform linear interpolation with clipping for out-of-bound data.
%
% Input arguments:
%   xBP: Break points on x-axis, monotonically increasing, size m x 1
%   V: Function values at break points, size m x 1
%   x: Interpolation points on x-axis, size q x 1
%
% Output arguments:
%   Vx: Interpolated points at x, size q x 1
%
% W. D. Widanage 04/06/2013 (Riffs!)

% Vectorize
xBP = xBP(:); 

% Get dimensions
m = length(xBP);
q = length(x);

% Initialise
Vx = zeros(q,1);

for qq = 1:q
    % Find indicies of nearest neighbours in x
    diffx = x(qq) - xBP;
    if (diffx < 0)
        ix0 = 1;
        ix1 = 1;
    elseif (diffx > 0)
        ix0 = m;
        ix1 = m;
    else
        diffx_0 = diffx;
        diffx_0(diffx < 0) = inf;
        [~,ix0] = min(diffx_0);
        diffx_1 = diffx;
        diffx_1(diffx > 0) = -inf;
        [~,ix1] = max(diffx_1);
    end
    
    x0 = xBP(ix0); 
    x1 = xBP(ix1);  

    if x0 == x1, xd = x0; else xd = (x(qq)-x0)/(x1-x0); end
       
    % Interplolate along x
    Vx(qq) = V(ix0)*(1-xd) + V(ix1)*xd;    

end
end



