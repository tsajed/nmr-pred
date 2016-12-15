% Converts axiality and rhombicity representation of the anisotro-
% pic part of a 3x3 interaction tensor into the corresponding mat-
% rix. Euler angles should be specified in radians. Syntax:
%
%                 M=axrh2mat(iso,ax,rh,alp,bet,gam)
%
% Parameters:
%
%        iso  - isotropic part of the interaction
%
%         ax  - interaction axiality, defined as zz-(xx+yy)/2
%               in terms of eigenvalues
%
%         rh  - interaction rhombicity, defined as (xx-yy) in
%               terms of eigenvalues
%
%        alp  - alpha Euler angle in radians
%
%        bet  - beta Euler angle in radians
%
%        gam  - gamma Euler angle in radians
%
% i.kuprov@soton.ac.uk

function M=axrh2mat(iso,ax,rh,alp,bet,gam)

% Check consistency
grumble(iso,ax,rh,alp,bet,gam);

% Compute eigenvalues
xx=(2*iso-2*ax+3*rh)/6;
yy=(2*iso-2*ax-3*rh)/6;
zz=(1*iso+2*ax)/3;

% Rotate the matrix
R=euler2dcm(alp,bet,gam);
M=R*diag([xx yy zz])*R';

end

% Consistency enforcement
function grumble(iso,ax,rh,alp,bet,gam)
if (~isnumeric(iso))||(~isreal(iso))||(~isscalar(iso))||...
   (~isnumeric(ax))||(~isreal(ax))||(~isscalar(ax))||...
   (~isnumeric(rh))||(~isreal(rh))||(~isscalar(rh))||...
   (~isnumeric(alp))||(~isreal(alp))||(~isscalar(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~isscalar(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~isscalar(gam))
    error('all inputs must be real scalars.');
end
end

% I'm working to improve my methods, and every hour
% I save is an hour added to my life.
%
% Ayn Rand, "Atlas Shrugged"

