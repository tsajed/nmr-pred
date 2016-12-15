% Pauli matrices for a spin of user-specified multiplicity. Syntax:
%
%                          sigma=pauli(mult)
%
% Where mult is an integer specifying the multiplicity and the following
% fields are returned in the output:
%
%    sigma.u - unit operator
%
%    sigma.p - raising operator
%
%    sigma.m - lowering operator
%
%    sigma.x - Pauli sigma_x matrix
%
%    sigma.y - Pauli sigma_y matrix
%
%    sigma.z - Pauli sigma_z matrix
%
% The matrices are normalized so as to fulfil the following commutation
% relations: 
%
%       [sigma.x,sigma.y]=1i*pauli.z
%       [sigma.y,sigma.z]=1i*pauli.x
%       [sigma.z,sigma.x]=1i*pauli.y
%
% with the raising and lowering operators defined as:
%
%        sigma.p=sigma.x+1i*sigma.y
%        sigma.m=sigma.x-1i*sigma.y
%
% i.kuprov@soton.ac.uk

function sigma=pauli(mult)

% Make sure the input is valid
if (~isnumeric(mult))||(~isscalar(mult))||(~isreal(mult))||...
   (mult<0)||(mod(mult,1)~=0)
    error('spin multiplicity must be a positive integer.');
end

% Get the Pauli matrices
if mult==2
    
    % Spin-half matrices are hard-coded for speed
    sigma.u=sparse([1 0; 0 1]);
    sigma.p=sparse([0 1; 0 0]);
    sigma.m=sparse([0 0; 1 0]);
    sigma.z=sparse([0.5 0; 0 -0.5]);
    sigma.x=0.5*(sigma.p+sigma.m);
    sigma.y=-0.5*1i*(sigma.p-sigma.m);
    
else
    
    % Everything else goes through the standard procedure
    spin=(mult-1)/2;
    prjs=((mult-1):-1:0)-spin;
    sigma.u=speye(mult,mult);
    sigma.p=spdiags(sqrt(spin*(spin+1)-prjs.*(prjs+1))',1,mult,mult);
    sigma.m=spdiags(sqrt(spin*(spin+1)-prjs.*(prjs-1))',-1,mult,mult);
    sigma.x=0.5*(sigma.p+sigma.m);
    sigma.y=-0.5*1i*(sigma.p-sigma.m);
    sigma.z=spdiags(prjs',0,mult,mult);
    
end
    
end

% Any refusal to recognize reality, for any reason whatever, has dis-
% astrous consequences. There are no evil thoughts except one -- the
% refusal to think. Don't ignore your own desires... Don't sacrifice
% them. Examine their cause.
%
% Ayn Rand, "Atlas Shrugged"

