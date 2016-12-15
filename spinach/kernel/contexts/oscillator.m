% Harmonic oscillator infrastructure. Syntax:
%
%           [H_oscl,X_oscl,xgrid]=oscillator(parameters)
%
% where:
%
%   parameters.frc_cnst    - force constant
%
%   parameters.dmd_mass    - particle mass
%
%   parameters.grv_cnst    - gravitational acceleration
%
%   parameters.n_points    - number of discretization points
%
%   parameters.rgn_size    - oscillator box dimension
%
% The function returns:
%
%   H_oscl   -   oscillator Hamiltonian
%
%   X_oscl   -   oscillator X operator
%
%   xgrid    - X coordinate grid
%
% i.kuprov@soton.ac.uk

function [H_oscl,X_oscl,xgrid]=oscillator(parameters)

% Second derivative operator
d2_dx2=((parameters.n_points-1)/parameters.rgn_size)^2*fdmat(parameters.n_points,5,2);

% Coordinate grid
xgrid=linspace(-parameters.rgn_size/2,parameters.rgn_size/2,parameters.n_points)';

% Coordinate operator
X_oscl=spdiags(xgrid,0,parameters.n_points,parameters.n_points);

% Boundary conditions
d2_dx2(:,1)=0; d2_dx2(:,end)=0; d2_dx2(1,:)=0; d2_dx2(end,:)=0;

% Hamiltonian
H_oscl=-(1/(2*parameters.dmd_mass))*d2_dx2+(parameters.frc_cnst/2)*X_oscl^2+...
       parameters.dmd_mass*parameters.grv_cnst*X_oscl;

end

% Just in terms of allocation of time resources, religion is
% not very efficient. There's a lot more I could be doing on
% a Sunday morning.
%
% Bill Gates

