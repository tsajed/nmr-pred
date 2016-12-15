% Converts the [RF_x, RF_y] representation of a pulse waveform and the
% derivatives of any function with respect to those RF values into the
% [RF_amplitude, RF_phase] representation and the derivatives of the
% function with respect to the amplitudes and phases. Syntax:
%
%        [A,phi,df_dA,df_dphi]=cartesian2polar(x,y,df_dx,df_dy)
%
% Parameters:
%
%    x       - vector of waveform amplitudes along X
%
%    y       - vector of waveform amplitudes along Y
%
%    df_dx   - optional vector of derivatives of a scalar function
%              with respect to the waveform amplitudes along X
%
%    df_dy   - optional vector of derivatives of a scalar function
%              with respect to the waveform amplitudes along Y
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
% Outputs:
%
%    A       - vector of waveform amplitudes
%
%    phi     - vector of waveform phases
%
%    df_dA   - vector of derivatives of the function with respect
%              to the waveform amplitudes.
%
%    df_dphi - vector of derivatives of the function with respect
%              to the waveform phases.
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
% i.kuprov@soton.ac.uk

function [A,phi,df_dA,df_dphi,d2f_dA2,d2f_dAdphi,d2f_dphi2]=cartesian2polar(x,y,df_dx,df_dy,d2f_dx2,d2f_dxdy,d2f_dy2)

% Check consistency
if nargin==2
    grumble(x,y);
elseif nargin==4
    grumble(x,y,df_dx,df_dy);
elseif nargin==7
    grumble(x,y,df_dx,df_dy,d2f_dx2,d2f_dxdy,d2f_dy2);
else
    error('incorrect number of arguments.');
end

% Transform coordinates
A=sqrt(x.^2+y.^2); phi=atan2(y,x);

% Transform derivatives
if (nargin>2)&&(nargout>2)
    df_dA=df_dx.*cos(phi)+df_dy.*sin(phi);
    df_dphi=-df_dx.*y+df_dy.*x;
end

% Transform second derivatives
if (nargin>4)&&(nargout>4)
    d2f_dphi2=  -A.*(cos(phi).*df_dx+sin(phi).*df_dy)+...
              A.*A.*(cos(phi).*cos(phi).*d2f_dy2-...
                   2*cos(phi).*sin(phi).*d2f_dxdy+...
                     sin(phi).*sin(phi).*d2f_dx2);
    d2f_dA2=cos(phi).*cos(phi).*d2f_dx2+...
          2*cos(phi).*sin(phi).*d2f_dxdy+...
            sin(phi).*sin(phi).*d2f_dy2;
    d2f_dAdphi=cos(phi).*df_dy-sin(phi).*df_dx+...
               A.*cos(phi).*sin(phi).*d2f_dy2+...
               A.*cos(phi).*cos(phi).*d2f_dxdy-...
               A.*sin(phi).*sin(phi).*d2f_dxdy-...
               A.*cos(phi).*sin(phi).*d2f_dx2;
end

end

% Consistency enforcement
function grumble(x,y,df_dx,df_dy,d2f_dx2,d2f_dxdy,d2f_dy2)
if nargin==2
    if (~isnumeric(x))||(~isreal(x))
        error('x parameter must be a vector of real numbers.');
    end
    if (~isnumeric(y))||(~isreal(y))
        error('y parameter must be a vector of real numbers.');
    end
    if ~all(size(x)==size(y))
        error('x and y vectors must have the same dimension.');
    end
elseif nargin==4
    if (~isnumeric(x))||(~isreal(x))
        error('x parameter must be a vector of real numbers.');
    end
    if (~isnumeric(y))||(~isreal(y))
        error('y parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dx))||(~isreal(df_dx))
        error('df_dx parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dy))||(~isreal(df_dy))
        error('df_dy parameter must be a vector of real numbers.');
    end
    if (~all(size(x)==size(y)))||(~all(size(y)==size(df_dx)))||(~all(size(df_dx)==size(df_dy)))
        error('all input vectors must have the same dimension.');
    end
elseif nargin==7
    if (~isnumeric(x))||(~isreal(x))
        error('x parameter must be a vector of real numbers.');
    end
    if (~isnumeric(y))||(~isreal(y))
        error('y parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dx))||(~isreal(df_dx))
        error('df_dx parameter must be a vector of real numbers.');
    end
    if (~isnumeric(df_dy))||(~isreal(df_dy))
        error('df_dy parameter must be a vector of real numbers.');
    end
    if (~isnumeric(d2f_dx2))||(~isreal(d2f_dx2))
        error('d2f_dx2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_dy2))||(~isreal(d2f_dy2))
        error('d2f_dy2 parameter must be a matrix of real numbers.');
    end
    if (~isnumeric(d2f_dxdy))||(~isreal(d2f_dxdy))
        error('ddf_dxdy parameter must be a matrix of real numbers.');
    end
    if (~all(size(x)==size(y)))||(~all(size(y)==size(df_dx)))||(~all(size(df_dx)==size(df_dy)))
        error('all input vectors must have the same dimension.');
    end
    if (size(d2f_dx2,2)~=length(df_dx))||(size(d2f_dx2,1)~=size(d2f_dx2,2))||all(size(d2f_dx2)~=size(d2f_dy2))||all(size(d2f_dy2)~=size(d2f_dxdy))
        error('all input matrices must have the same, square dimensions.');
    end
end
end

% Beware of geeks bearing formulas.
%
% Warren Buffett 

