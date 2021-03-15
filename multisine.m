function varargout=multisine(varargin)
% Johan Paduart. Last Updated: 18/01/2006
% Creates a random phase multisine. 
% 
% Use:  u = multisine(N, kind, M, R {, Nblock});
% N         : Number of points
% kind      : 'Full','Odd','SpecialOdd', or 'RandomOdd'
% M         : Last excited line
% R         : Number of Realisations
% Nblock    : Block Length for Random Odd Multisine
% Default values (without arguments): N=1024; 'kind'='Full'; M=90%
% 
% Output arguments:
% u = multisine(...)
% [u,lines] = multisine(...)
% [u,lines,non_exc] = multisine(...)
% [u,lines,non_exc_odd,non_exc_even] = multisine(...)
% 
% Another way to use multisine is:
% u = multisine(N,lines);
% e.g. when you want to create another realisation, with exactly the same 
% excited lines.

% Default values:
N=1024;kind='full';M=floor(0.9*N/2);R=1;Nblock=4;

switch nargin
    case 0 
        % No input arguments => use default values
    case 1 % Only one argument => number of points
        N=varargin{1};
        M=floor(0.9*N/2);
    case 2 % Two arguments => number of points and kind OR number and lines
        [N,temp]=deal(varargin{:});
        M=floor(0.9*N/2);
        if isa(temp,'char')
            kind=temp;
        else
            kind='dummy'; lines=temp;
        end
    case 3 % Nbr of points, kind and bandwidth
        [N,temp1,temp2]=deal(varargin{:});
        if isa(temp1,'char')
            kind=temp1;
            M=temp2;
        else
            kind='dummy'; lines=temp1; R=temp2;
        end
    case 4 % Only when Random Odd Multisine
        [N,kind,M,R]=deal(varargin{:});
    case 5
        [N,kind,M,R,Nblock]=deal(varargin{:});
    otherwise
        error('Too many parameters')
end

switch lower(kind)
    case 'full'
        lines=2:M;
    case 'odd'
        lines=2:2:M;
    case {'specialodd', 'special odd','special'}
        lines= sort([2:8:M,4:8:M]);
    case {'randomodd','random odd','random'}
        lines=2:2:M;
        Nblocks=floor(length(lines)/Nblock);
        indices=1:Nblock:Nblocks*Nblock; % start indices of blocks
        indices=indices+floor(Nblock*rand(1,Nblocks)); % plus random nbr
        lines(indices)=[]; % eliminate selected detection lines
    case 'dummy'
        % do nothing
    otherwise
        error('Unknown type');
end
y=zeros(N,R);
y(lines,:)=1;
y=y.*exp(j*2*pi*rand(size(y)));
y=2*real(ifft(y));
switch nargout
    case 1
        varargout={y};
    case 2
        varargout={y,lines};
    case 3
        non_lines=1:N/2;
        non_lines(lines)=[];
        varargout={y,lines,non_lines};
    case 4
        non_lines=1:N/2;
        non_lines(lines)=[];
        non_odd=non_lines(logical(mod(non_lines-1,2)));
        non_even=non_lines(~mod(non_lines-1,2));
        varargout={y,lines,non_odd,non_even};
end
