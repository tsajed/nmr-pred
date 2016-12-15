% Converts the C_q and eta_q quadrupolar interaction specification con-
% vention into a 3x3 interaction matrix in Hz. Usage:
%
%                    Q=eeqq2nqi(C_q,eta_q,I,euler_angles)
%
% where C_q is the quadrupolar coupling constant (Hz), eta_q is the quad- 
% rupolar tensor asymmetry parameter, I is the spin quantum number and
% the last parameter is a vector of three Euler angles in radians, giving
% the orientation of the principal axis frame relative to the lab frame.
%
% i.kuprov@soton.ac.uk

function Q=eeqq2nqi(C_q,eta_q,I,eulers)

% Check consistency
grumble(C_q,eta_q,I,eulers);

% Get the eigenvalues
XX=-C_q*(1-eta_q)/(4*I*(2*I-1));
YY=-C_q*(1+eta_q)/(4*I*(2*I-1));
ZZ=+C_q/(2*I*(2*I-1));

% Get the rotation matrix
R=euler2dcm(eulers);

% Get the quadrupole tensor matrix
Q=R*diag([XX YY ZZ])*R';

end

% Consistency enforcement
function grumble(C_q,eta_q,I,eulers)
if (~isnumeric(C_q))||(~isnumeric(eta_q))||(~isnumeric(I))||(~isnumeric(eulers))
    error('all inputs must be numeric.');
end
if (~isreal(C_q))||(~isreal(eta_q))||(~isreal(I))||(~isreal(eulers))
    error('all inputs must be real.');
end
if numel(eulers)~=3
    error('eulers vector must have three elements.');
end
if (numel(C_q)~=1)||(numel(eta_q)~=1)
    error('C_a and eta_q arguments must have a single element.');
end
if (numel(I)~=1)||(I<1)||(mod(2*I+1,1)~=0)
    error('I must be an integer or half-integer greater or equal to 1.');
end
end

% A true man does what he will, not what he must. 
%
% George R.R. Martin, "Game of Thrones"

