% Computes the dihedral angle between vectors specified
% by the four sets of atomic coordinates. The angle is
% returned in degrees. The atoms are assumed to be bon-
% ded as A-B-C-D. Syntax:
%
%                 phi=dihedral(A,B,C,D)
%
% i.kuprov@soton.ac.uk

function phi=dihedral(A,B,C,D)

% Check consistency
grumble(A,B,C,D);

% Do the math
b1=B-A; b2=C-B; b3=D-C; b1=b1/norm(b1); b2=b2/norm(b2); b3=b3/norm(b3);
phi=180*atan2(dot(norm(b2)*b1,cross(b2,b3)),dot(cross(b1,b2),cross(b2,b3)))/pi;

end

% Consistency enforcement
function grumble(A,B,C,D)
if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||(~isnumeric(D))||...
   (~isreal(A))||(~isreal(B))||(~isreal(C))||(~isreal(D))||...
   (numel(A)~=3)||(numel(B)~=3)||(numel(C)~=3)||(numel(D)~=3)
    error('the arguments must be 3-element vectors of real numbers.');
end
end

% IK's bitcoin address is 18p5ttXFwyqbtiLAhUCDWk2vo61DVUZS21. He
% probably deserves a pint of beer, right? ;)

