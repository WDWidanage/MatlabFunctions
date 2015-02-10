function [sol, fsol, counter] = BrentRootMethod(fcnHandle,a,b,errorTol)
%
% Finds the root of a function f(x) = 0 via the Brent's root method.
% Root should lie between a and b and f(a) and f(b) have opposite sign.
%
% Input arguments
%   fcnHandle : Function handle of function to solve, fcnHandle = @fcn(...)
%   a         : Initial lower limit, size 1 x 1
%   b         : Initial upper limit, size 1 x 1
%   errorTol  : Error tolerance of root, positive real value, size 1 x 1
% Output argument
%   sol: root, size 1 x 1
%   fsol: function evaluation at sol, size 1 x 1
%   counter: Number of iterations taken to reach sol, size 1 x 1
%
% W.D. Widanage 19/05/2013 (Spicy)


fa = fcnHandle(a);
fb = fcnHandle(b);

%   Initialise counter
counter = 1;

% If function does not have a sign change over intitial interval return the
% value closest to zero
if (fa*fb) >= 0
    if (abs(fa) < (fb))
        sol = a;
        fsol =fa;
        return;
    else
        sol= b;
        fsol =fb;
        return;
    end
end

% If |f(a)| < |f(b)| swap a, b and fa and fb
if abs(fa) < abs(fb)
    temp = a; a = b; b= temp;
    temp = fa; fa = fb; fb = temp;
end

c = a;
fc = fa;
mflag = true;


while (fb~=0) && (abs(a-b)) > errorTol
    if fa~= fc && fb ~= fc
        % Inverse quadratic interpolation
        s = a*fb*fc/(fa-fb)/(fa-fc)+b*fa*fc/(fb-fa)/(fb-fc)+c*fa*fb/(fc-fa)/(fc-fb);
    else
        % Secant Rule
        s = b-fb*(b-a)/(fb-fa);
    end
    
    temp = (3*a+b)/4;
    if ((s < temp) && (s < b)) || ((s > temp) && (s > b))...        % Condition 1: s is not between (3a + b)/4 and b
            || (mflag && (abs(s-b) >= (abs(b-c)/2)))...             % Condition 2: mflag is set and |s?b| ? |b?c| / 2
            || ((mflag == false) && (abs(s-b) >= (abs(c-d)/2)))...  % Condition 3: mflag is cleared and |s?b| ? |c?d| / 2
            || (mflag && (abs(b-c) < errorTol))...                  % Condition 4: mflag is set and |b?c| < errorTol
            || ((mflag == false) && (abs(c-d) < errorTol)),         % Condition 5: mflag is cleared and |c?d| < errorTol
        
        s = (a+b)/2;
        mflag = true;
    else
        mflag = false;
    end
    fs = fcnHandle(s);
    d = c;
    c = b;
    fc = fb;
    
    % Preserve ordinates that result in a sign change
    if (fa*fs) < 0
        b = s; fb = fs;
    else
        a = s; fa = fs;
    end
    
    % If |f(a)| < |f(b)| swap a, b and fa and fb
    if abs(fa) < abs(fb)
        temp = a; a = b; b= temp;
        temp = fa; fa = fb; fb = temp;
    end
    
    % Return b if counter is > 1000
    if counter > 1000
        sol= b;
        fsol = fb;
        return;
    end
    counter = counter +1;
end

sol = b;
fsol = fb;

