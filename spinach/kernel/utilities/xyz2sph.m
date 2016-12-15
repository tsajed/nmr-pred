% Converts Cartesian coordinates [x y z] into spherical coordinates
% according to the ISO convention.
%
% e.suturina@soton.ac.uk

function [r, theta, phi] = xyz2sph(x, y, z)

% Radius 0 <= r < Inf
r=sqrt(x.^2+y.^2+z.^2);

% Inclination 0 <= theta <= pi
theta=acos(z./r); 

% Azimuth 0 <= phi < 2*pi
phi=atan2(y,x);    

end

% Q: How would a theoretical physicist milk a cow?
% A: As a first approximation, consider a spherical cow in vacuum...

 