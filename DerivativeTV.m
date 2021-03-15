function [du,iterUpdate,iterPCG,flag] = DerivativeTV(u,x,alpha,iterMax,iterPCGMax,dia)
%
% Computes the first order derivative of the data signal u. If using a finite
% difference method it will amplify all the small changes (noise) in the u
% vector. Solution is to regularise the differentiation process using
% Total Variation regularisation.
%
% Mandotory input arguments:
%   u: Input signal to be differentiated, size N x 1
%
% Optional arguments
%   x:          Abscissa values at which u occurs (normally x is the time points of u), default x = [0 : 1/N : 1-1/N], size N x 1
%   alpha:      Regularisation scalar factor, default 1E-3, size 1 x 1
%   iterMax:    Number for iterations for updating the solution, default 100, size 1 x 1
%   iterPCGMax: Maximum number for iterations for the Preconditioned Congugate Gradient algorithm, default 100, size 1 x 1
%   dia:        Set dia to 1 for diagnosis and to 0 for no diagnosis
%
% Output arguments
%   du:         Regularised derivative of u, size N x 1
%   iterUpdate: Number of iterations of the soulution update, size 1 x 1
%   iterPCG:    Number of iterations for the PCG algorithm, size 1 x 1
%   flag:       Gives more info on PCG convergence details. See pcg help for flag
%               details.
%
% Adapted from "Numerical differentiation of noisy, nonsmooth data" 
% Applied Mathematics, Vol. 2011, Article ID 164564, 2011. paper and code.
% W. D. Widanage 04/09/2014 (They see me roooolling!)

u = u(:);               % Vectorise input data points
N  = length(u);         % Length of data signal

if nargin < 2 || isempty(x)     % Calculate horizontal spacing of the u data points
    dx = 1/N*ones(N,1);
    normFactor = 1;
else
    x = x(:);               % Vectorise abscissa data points
    % Normalise x vector
    normFactor = max(x-x(1));
    xNorm = x/normFactor;
    dx = diff(xNorm);
    dx = [dx(1); dx];
end
if nargin < 3 || isempty(alpha)
    alpha = 1E-3;
end
if nargin < 4 || isempty(iterMax)
    iterMax = 100;
end
if nargin < 5 || isempty(iterPCGMax)
    iterPCGMax = 100;  % This is the maximum number of iterations pcg can go up to if flag is still 1
end
if nargin < 6 || isempty(dia)
    dia = 0;
end

u0 = u - u(1); % Deduct first value from signal to ensure zero initial condition

% Initial solution: Proportional to finite difference derivative
m = [0; diff(u0)];

J = @(v) cumsum(v);                              % Represent the Jacobian matrix times m, which is simply the integration of the signal, as a function
JT = @(w) flipud(cumsum(flipud(w(:))));          % Represent Jacobian transpose times Jacobian times m as a function

% JT times data, required for gradient of function
JTu = JT(u0);

% Construct first-order finite difference matrix L associated with the regularisation term
e = ones(N,1);
L = spdiags(1./dx,0,N,N)*spdiags([-e,e],[-1,0],N,N);
L(1,1) = 0;
ep = 1e-8;    % Threshold value to use when dividing by zero
deltaM = 100; % Initialise solution update to a "large" value to start the loop
iterUpdate = 1;

% Approximate derivative by backward difference, derivative of first point
% should be ignored
finiteDiffU = spdiags(1./dx,0,N,N)*m;

% Start solution update
while norm(deltaM) > 1E-4 && iterUpdate <= iterMax
    
    % Update W matrix
    wi = 1./sqrt((L*m).^2 + ep);
    W = spdiags(wi,0,N,N);
    
    % Gradient of function with regularisation
    delF = JT(J(m)) - JTu + alpha*L'*W*L*m;
    
    % Hessian function
    H = @(x)JT(J(x)) + alpha*L'*W*L*x;
    
    % Preconditioning matrix
    c = cumsum(N:-1:1).';
    M = alpha*L'*W*L + spdiags(c(end:-1:1),0,N,N);
    opts.type = 'ict'; opts.droptol = 1e-2; %opts.shape = 'upper';
    R = ichol(M, opts);
    
    % Call pcg method
    tol = 1.0e-4;
    [deltaM,flag,relRes,iterPCG,resVec] = pcg(H,-delF,tol,iterPCGMax,R,R');
    if flag ~=0 && iterUpdate == iterMax
        warning(['PCG did not converged to the desired tolerance tol: ', num2str(tol),' within iterPCGMax: ',num2str(iterPCGMax),' iterations. Increase iterPCGMax argument'])
    end
    
    m = m + deltaM;
    
    % If diagnosis argument true plot in loop
    if dia
        if numel(dx)==1
            relTime = [0:N-1]*dx*normFactor;
        else
            relTime = cumsum(dx);
            relTime = (relTime - relTime(1))*normFactor;             % Relative horizontal vector now starts at zero
        end
        
        % Scale to get derivative solution
        du = spdiags(1./dx,0,N,N)*m;
        
        figure(100);
        subplot(3,1,1)
        plot(relTime,u);
        ylabel('Signal to differentiate u')
        subplot(3,1,2)
        plot(relTime,finiteDiffU);
        ylabel('Simple finite difference, dudx');
        subplot(3,1,3)
        plot(relTime,du); grid on;
        xlabel('x'); ylabel('TV regularised solution, dudx');
        title(['Update iteration no:', num2str(iterUpdate), ' Solution update norm: ',  num2str(chop(norm(deltaM),3)),' Gradient norm: ', num2str(chop(norm(delF),3))])
    end
    
    iterUpdate = iterUpdate+1;
end
% Scale to get derivative estimate
du = spdiags(1./dx,0,N,N)*m/normFactor;



