function Vxyz= LinearInter3D(xBP,yBP,zBP,V,x,y,z)

% Perform trilinear interpolation with clipping for out-of-bound data.
%
% Input arguments:
%   xBP: Break points on x-axis, monotonically increasing, size m x 1
%   yBP: Break points on y - axis, monotonically increasing, size n x 1
%   zBP: Break points on z-axis, monotonically increasing, size p x 1
%   V: Function values at break points, size m x n x p
%   x: Interpolation points on x-axis, size q x 1
%   y: Interpolation points on y-axis, size q x 1
%   z: Interpolation points on z-axis, size q x 1
%
% Output arguments:
%   Vxyz: Interpolated points at x,y and z, size q x 1
%
% W. D. Widanage 03/06/2013 (hungry!)

% Vectorize
xBP = xBP(:); yBP = yBP(:); zBP = zBP(:);

% Get dimensions
m = length(xBP);
n = length(yBP);
p = length(zBP);
q = length(x);

% Initialise
Vxyz = zeros(q,1);

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
    % Find indicies of nearest neighbours in z
    diffz = z(qq) - zBP;
    if (diffz < 0)
        iz0 = 1;
        iz1 = 1;
    elseif (diffz > 0)
        iz0 = p;
        iz1 = p;
    else
        diffz_0 = diffz;
        diffz_0(diffz < 0) = inf;
        [~,iz0] = min(diffz_0);
        diffz_1 = diffz;
        diffz_1(diffz > 0) = -inf;
        [~,iz1] = max(diffz_1);
    end
    
    x0 = xBP(ix0); y0 = yBP(iy0); z0 = zBP(iz0);
    x1 = xBP(ix1); y1 = yBP(iy1); z1 = zBP(iz1);

    if x0 == x1, xd = x0; else xd = (x(qq)-x0)/(x1-x0); end
    if y0 == y1, yd = y0; else yd = (y(qq)-y0)/(y1-y0); end
    if z0 == z1, zd = z0; else zd = (z(qq)-z0)/(z1-z0); end
       
    % Interplolate along x
    V00 = V(ix0,iy0,iz0)*(1-xd) + V(ix1,iy0,iz0)*xd;
    V10 = V(ix0,iy1,iz0)*(1-xd) + V(ix1,iy1,iz0)*xd;
    V01 = V(ix0,iy0,iz1)*(1-xd) + V(ix1,iy0,iz1)*xd;
    V11 = V(ix0,iy1,iz1)*(1-xd) + V(ix1,iy1,iz1)*xd;
    
    % Interpolate along y
    V0 = V00*(1-yd)+V10*yd;
    V1 = V01*(1-yd)+V11*yd;
    
    % Interpolate along z
    Vxyz(qq) = V0*(1-zd)+V1*zd;
end
end

