function Vxy= LinearInter2D(xBP,yBP,V,x,y)

% Perform bilinear interpolation with clipping for out-of-bound data.
%
% Input arguments:
%   xBP: Break points on x-axis, monotonically increasing, size m x 1
%   yBP: Break points on y - axis, monotonically increasing, size n x 1
%   V: Function values at break points, size m x n
%   x: Interpolation points on x-axis, size q x 1
%   y: Interpolation points on y-axis, size q x 1
%
% Output arguments:
%   Vxy: Interpolated points at x and y, size q x 1
%
% W. D. Widanage 04/06/2013 (Riffs!)

% Vectorize
xBP = xBP(:); yBP = yBP(:);

% Get dimensions
m = length(xBP);
n = length(yBP);
q = length(x);

% Initialise
Vxy = zeros(q,1);

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
    % Find indicies of nearest neighbours in y
    diffy = y(qq) - yBP;
    if (diffy < 0)
        iy0 = 1;
        iy1 = 1;
    elseif (diffy > 0)
        iy0 = n;
        iy1 = n;
    else
        diffy_0 = diffy;
        diffy_0(diffy < 0) = inf;
        [~,iy0] = min(diffy_0);
        diffy_1 = diffy;
        diffy_1(diffy > 0) = -inf;
        [~,iy1] = max(diffy_1);
    end
    
    x0 = xBP(ix0); y0 = yBP(iy0); 
    x1 = xBP(ix1); y1 = yBP(iy1); 

    if x0 == x1, xd = x0; else xd = (x(qq)-x0)/(x1-x0); end
    if y0 == y1, yd = y0; else yd = (y(qq)-y0)/(y1-y0); end
       
    % Interplolate along x
    V0 = V(ix0,iy0)*(1-xd) + V(ix1,iy0)*xd;
    V1 = V(ix0,iy1)*(1-xd) + V(ix1,iy1)*xd;
    
    % Interpolate along z
    Vxy(qq) = V0*(1-yd)+V1*yd;
end
end



