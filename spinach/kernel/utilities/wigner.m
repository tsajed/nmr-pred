% Computes Wigner matrices of user-specified ranks. Syntax:
%
%                  D=wigner(l,alp,bet,gam)
%
% where alp, bet and gam are Euler angles in radians. Rows and columns
% of the resulting Wigner matrix are sorted by descending ranks, e.g.:
%    
%                         [D( 2,2)  ...  D( 2,-2)
%                            ...    ...    ...  
%                          D(-2,2)  ...  D(-2,-2)]
%
% The resulting Wigner matrix is to be used as v=W*v, where v is a co-
% lumn vector of irreducible spherical tensor coefficients, listed ver-
% tically in the order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).
%
% i.kuprov@soton.ac.uk

function D=wigner(l,alp,bet,gam)

% Check consistency
grumble(l,alp,bet,gam);

% Get Pauli matrices
L=pauli(2*l+1);

% Compute Wigner matrix (Brink & Satchler, Eq 2.13)
D=expm(-1i*L.z*alp)*expm(-1i*L.y*bet)*expm(-1i*L.z*gam);

end

% Consistency enforcement
function grumble(l,alp,bet,gam)
if (~isnumeric(l))||(~isreal(l))||(~isscalar(l))||(mod(l,1)~=0)||(l<0)
    error('l must be a non-negative real integer.');
end
if (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('alp, bet and gam must be real scalars.');
end
end

% Every stink that fights the ventilator thinks it
% is Don Quixote.
%
% Stanistaw Lec

