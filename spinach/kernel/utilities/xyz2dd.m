% Converts coordinate specification of the dipolar interaction
% into the dipolar interaction constant (angular frequency un-
% its), three Euler angles in radians and the dipolar interac-
% tion matrix (angular frequency units). Syntax:
%
%      [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)
%
% where r1 and r2 are 3-element vectors of spin coordinates in
% Angstroms and isotope(1,2) are isotope strings, e.g. '13C'.
%
% N.B. Euler angles are not uniquely defined for the orientati-
%      on of axial interactions (gamma angle can be anything).
%
% i.kuprov@soton.ac.uk

function [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)

% Fundamental constants
hbar=1.054571628e-34; mu0=4*pi*1e-7;

% Get the distance
distance=norm(r2-r1,2);
        
% Get the ort
ort=(r2-r1)/distance;
        
% Get the dipolar interaction constant
d=spin(isotope1)*spin(isotope2)*hbar*mu0/(4*pi*(distance*1e-10)^3);
        
% Get the Euler angles
[alp,bet,~]=cart2sph(ort(1),ort(2),ort(3)); bet=pi/2-bet; gam=0;

% Get the dipolar coupling matrix and symmetrise it
m=d*[1-3*ort(1)*ort(1)   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
      -3*ort(2)*ort(1)  1-3*ort(2)*ort(2)   -3*ort(2)*ort(3);
      -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)*ort(3)];
m=(m+m')/2;

end

% This principle is not a theorem, but a physical proposition, 
% that is, a vaguely stated and, strictly speaking, false as-
% sertion. Such assertions often happen to be fruitful sourc-
% es for mathematical theorems.
%
% Vladimir Arnold

