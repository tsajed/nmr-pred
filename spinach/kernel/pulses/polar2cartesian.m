% Converts [RF_amplitude, RF_phase] representation of a pulse waveform and
% the derivatives of any function with respect to those amplitudes and pha-
% ses into the [RF_x, RF_y] representation and the derivatives of the func-
% tion with respect to those X and Y RF values. Syntax:
%
%         [x,y,df_dx,df_dy]=polar2cartesian(A,phi,df_dA,df_dphi)
%
% Parameters:
%
%    A       - vector of waveform amplitudes
%
%    phi     - vector of waveform phases
%
%    df_dA   - optional vector of derivatives of some scalar function  
%              with respect to the waveform amplitudes.
%
%    df_dphi - optional vector of derivatives of some scalar function
%              with respect to the waveform phases.
%
%    d2f_dA2 - matrix of second derivatives of the function with respect
%              to the waveform amplitudes.
%
%  d2f_dAdphi- matrix of second derivatives of the function with respect
%              to the waveform amplitudes and phases.
%
%   d2f_dphi2- matrix of second derivatives of the function with respect
%              to the waveform phases.
%
% Outputs:
%
%    x       - vector of waveform amplitudes along X
%
%    y       - vector of waveform amplitudes along Y
%
%    df_dx   - vector of derivatives of the function with respect to 
%              the waveform amplitudes along X
%
%    df_dy   - vector of derivatives of the function with respect to 
%              the waveform amplitudes along Y
%
%    d2f_dx2 - optional matrix of second derivatives of a scalar function
%              with respect to the waveform amplitudes along X
%
%    d2f_dxdy- optional matrix of second derivatives of a scalar function
%              with respect to the waveform amplitudes along X and Y
%
%    d2f_dy2 - optional matrix of second derivatives of a scalar function
%              with respect to the waveform amplitudes along Y
%
% i.kuprov@soton.ac.uk

function [x,y,df_dx,df_dy,d2f_dx2,d2f_dxdy,d2f_dy2]=polar2cartesian(A,phi,df_dA,df_dphi,d2f_dA2,d2f_dAdphi,d2f_dphi2)

% Check consistency
if nargin==2
    grumble(A,phi);
elseif nargin==4
    grumble(A,phi,df_dA,df_dphi);
elseif nargin==7
    grumble(A,phi,df_dA,df_dphi,d2f_dA2,d2f_dAdphi,d2f_dphi2);
else
    error('incorrect number of arguments.');
end

% Transform coordinates
x=A.*cos(phi); y=A.*sin(phi);

% Transform derivatives
if (nargin>2)&&(nargout>2)
    df_dx=df_dA.*cos(phi)-df_dphi.*sin(phi)./A;
    df_dy=df_dA.*sin(phi)+df_dphi.*cos(phi)./A;
end

% Transform second derivatives
if (nargin>4)&&(nargout>4)
    A4=A.^4; A3=A.^3; A2=A.^2;
    d2f_dx2=+2*x.*y.*df_dphi./A4+y.*y.*d2f_dphi2./A4+df_dA./A-x.*x.*df_dA./A3-2*x.*y.*d2f_dAdphi./A3+x.*x.*d2f_dA2./A2;
    d2f_dy2=-2*x.*y.*df_dphi./A4+x.*x.*d2f_dphi2./A4+df_dA./A-y.*y.*df_dA./A3+2*x.*y.*d2f_dAdphi./A3+y.*y.*d2f_dA2./A2;
    d2f_dxdy=-df_dphi./A2+2*y.*y.*df_dphi./A4-x.*y.*d2f_dphi2./A4-x.*y.*df_dA./A3+x.*x.*d2f_dAdphi./A3-y.*y.*d2f_dAdphi./A3+x.*y.*d2f_dA2./A2;
end

end

% Consistency enforcement
function grumble(A,phi,df_dA,df_dphi,d2f_dA2,d2f_dAdphi,d2f_dphi2)
if nargin==2
    if (~isnumeric(A))||(~isreal(A))||(~all(A>=0))
        error('amplitude parameter must be a vector of non-negative real numbers.');
    end
    if (~isnumeric(phi))||(~isreal(phi))
        error('phase parameter must be a vector of real numbers.');
    end
    if ~all(size(A)==size(phi))
        error('amplitude and phase vectors must have the same dimension.');
    end
elseif nargin==4
    if (~isnumeric(A))||(~isreal(A))||(~all(A>=0))
        error('amplitude parameter must be a vector of non-negative real numbers.');
    end
    if (~isnumeric(phi))||(~isreal(phi))
        error('phase parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dA))||(~isreal(df_dA))
        error('df_dA parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dphi))||(~isreal(df_dphi))
        error('df_dphi parameter must be a vector of real numbers.');
    end
    if (~all(size(A)==size(phi)))||(~all(size(phi)==size(df_dA)))||(~all(size(df_dA)==size(df_dphi)))
        error('all input vectors must have the same dimension.');
    end
elseif nargin==7
    if (~isnumeric(A))||(~isreal(A))||(~all(A>=0))
        error('amplitude parameter must be a vector of non-negative real numbers.');
    end
    if (~isnumeric(phi))||(~isreal(phi))
        error('phase parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dA))||(~isreal(df_dA))
        error('df_dA parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dphi))||(~isreal(df_dphi))
        error('df_dphi parameter must be a vector of real numbers.');
    end
    if (~all(size(A)==size(phi)))||(~all(size(phi)==size(df_dA)))||(~all(size(df_dA)==size(df_dphi)))
        error('all input vectors must have the same dimension.');
    end
    if (~isnumeric(d2f_dA2))||(~isreal(d2f_dA2))
        error('d2f_dA2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_dphi2))||(~isreal(d2f_dphi2))
        error('d2f_dphi2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_dAdphi))||(~isreal(d2f_dAdphi))
        error('ddf_dAdphi parameter must be a matrix of real numbers.');
    end
    if (size(d2f_dA2,2)~=length(df_dA))||(size(d2f_dA2,1)~=size(d2f_dA2,2))||all(size(d2f_dA2)~=size(d2f_dphi2))||all(size(d2f_dphi2)~=size(d2f_dAdphi))
        error('all input matrices must have the same, square dimensions.');
    end
end
end

% A public opinion poll is no substitute for thought.
%
% Warren Buffett

