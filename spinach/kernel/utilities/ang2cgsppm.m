% Converts magnetic susceptibility from the Angstrom^3 units 
% required by Spinach pseudocontact shift functionality into
% the cgs-ppm (aka cm^3/mol) units quoted by quantum chemist-
% ry packages into. Syntax:
%
%                     cgsppm=ang2cgsppm(ang)
%
% Arrays of any dimension are supported.
%
% i.kuprov@soton.ac.uk

function cgsppm=ang2cgsppm(ang)

% Check consistency
grumble(ang);

% Do the calculation
cgsppm=6.02214129e23*ang/(4*pi*1e18);

end

% Consistency enforcement
function grumble(ang)
if ~isnumeric(ang)
    error('input must be numeric.');
end
end

% No artist tolerates reality.
%
% Friedrich Nietzsche

