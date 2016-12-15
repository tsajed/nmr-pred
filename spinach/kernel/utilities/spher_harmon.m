% Computes spherical harmonics. Syntax:
%
%                      Y=spher_harmon(l,m,theta,phi)
%
% Parameters:
%
%     l        - L quantum number
%
%     m        - M quantum number
%
%     theta    - a vector of theta angles in radians
%
%     phi      - a vector of phi angles in radians
%
% i.kuprov@soton.ac.uk

function Y=spher_harmon(l,m,theta,phi)

% Check consistency
grumble(l,m,theta,phi);

% Get Schmidt-normalized Legendres
S=legendre(l,cos(theta),'sch');
S=S(abs(m)+1,:); S=S(:);

% Make spherical harmonics
if m==0
    Y=sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi);
else
    Y=sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi)/sqrt(2);
end

% Flip the sign if needed
if (m>0)&&(mod(m,2)==1), Y=-Y; end

end

% Consistency enforcement
function grumble(l,m,theta,phi)
if (~isnumeric(l))||(~isreal(l))||(~isscalar(l))||(mod(l,1)~=0)||(l<0)
    error('l must be a non-negative real integer.');
end
if (~isnumeric(m))||(~isreal(m))||(~isscalar(m))||(mod(m,1)~=0)||(m<-l)||(m>l)
    error('m must be a real integer from [-l,l] interval.');
end
if (~isnumeric(theta))||(~isreal(theta))||(size(theta,2)~=1)
    error('theta must be a real scalar or column vector.');
end
if (~isnumeric(phi))||(~isreal(phi))||(size(phi,2)~=1)
    error('phi must be a real scalar or column vector.');
end
end

% Don't pay any attention to what they write about 
% you. Just measure it in inches.
%
% Andy Warhol

