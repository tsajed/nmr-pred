% Generates spatially distributed dynamics operators. Syntax:
%
%  [H,R,K,D,Gx,Gy,Gz,Fx,Fy,Fz]=imaging(H,Ph,K,Z,parameters)
%
% The following outputs are returned: 
%
%    H        - spin Hamiltonian in the Fokker-Planck space
%
%    R        - spin relaxation superoperator in the Fokker-
%               Planck space
%
%    K        - chemical kinetics superoperator in the Fokk-
%               er-Planck space
%
%    Gx,Gy,Gz - gradient operators (normalized for 1.0 T/m)
%
%    D        - diffusion operator (absorptive boundary con-
%               conditions, normalized to 1.0 m^2/s)
%
% The following parameters must be provided:
%
%    H        - spin Hamiltonian commutation superoperator 
%               in Liouville space
%
%    Ph       - relaxation superoperators and their phantoms
%               as {{R1,Ph1},{R2,Ph2},...}
%
%    K        - chemical kinetics superoperator in Liouvil-
%               le space
%
%    Z        - normalized Zeeman Hamiltonian commutation
%               superoperator in Liouville space(to be ob-
%               tained from carrier.m and divided by the 
%               magnet field)
%
%   parameters.dims    - dimensions of the 3D box, meters
%
%   parameters.npts    - number of points in each dimension
%                        of the 3D box
%
%   parameters.cond    - {'fourier'} for Furier derivative operators;
%                        {'period',n} for n-point finite-difference
%                        operators with periodic boundary condition
%
% Note: gradients are assumed to be linear and centered on the
%       middle of the sample.
%
% Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
%       responds to a column-wise vectorization of a 3D array
%       with dimensions ordered as [X Y Z].
%
% i.kuprov@soton.ac.uk
% a.j.allami@soton.ac.uk

function [H,R,K,D,Gx,Gy,Gz,Fx,Fy,Fz]=imaging(H,Ph,K,Z,parameters)

% Adapt to the dimensionality
switch numel(parameters.npts)
    
    case 1
        
        % Call the hydrodynamics infrastructure
        Fx=hydrodynamics(parameters); Fy=[]; Fz=[];
        
        % Get the diffusion operator
        D=-1i*Fx*Fx;
        
        % Kron into Fokker-Planck space
        Fx=kron(Fx,speye(size(Z)));
        D=kron(D,speye(size(Z)));
        
        % Generate normalized gradient operators
        Gx=linspace(-0.5,0.5,parameters.npts(1));
        Gx=spdiags(Gx',0,parameters.npts(1),parameters.npts(1));
        Gx=parameters.dims(1)*kron(Gx,Z); Gy=[]; Gz=[];
        
    case 2
        
        % Call the hydrodynamics infrastructure
        [Fx,Fy]=hydrodynamics(parameters); Fz=[];
        
        % Get the diffusion operator
        D=-1i*(Fx*Fx+Fy*Fy);
        
        % Kron into Fokker-Planck space
        Fx=kron(Fx,speye(size(Z)));
        Fy=kron(Fy,speye(size(Z)));
        D=kron(D,speye(size(Z)));
        
        % Generate normalized gradient operators
        IdX=speye(parameters.npts(1)); 
        IdY=speye(parameters.npts(2));
        Gx=linspace(-0.5,0.5,parameters.npts(1));
        Gy=linspace(-0.5,0.5,parameters.npts(2));
        Gx=spdiags(Gx',0,parameters.npts(1),parameters.npts(1));
        Gy=spdiags(Gy',0,parameters.npts(2),parameters.npts(2));
        Gx=parameters.dims(1)*kron(IdY,kron(Gx,Z)); 
        Gy=parameters.dims(2)*kron(Gy,kron(IdX,Z)); Gz=[];
        
    case 3
        
        % Call the hydrodynamics infrastructure
        [Fx,Fy,Fz]=hydrodynamics(parameters);
        
        % Get the diffusion operator
        D=-1i*(Fx*Fx+Fy*Fy+Fz*Fz);
        
        % Kron into Fokker-Planck space
        Fx=kron(Fx,speye(size(Z)));
        Fy=kron(Fy,speye(size(Z)));
        Fz=kron(Fz,speye(size(Z)));
        D=kron(D,speye(size(Z)));
        
        % Generate normalized gradient operators
        IdX=speye(parameters.npts(1));
        IdY=speye(parameters.npts(2));
        IdZ=speye(parameters.npts(3));
        Gx=linspace(-0.5,0.5,parameters.npts(1));
        Gy=linspace(-0.5,0.5,parameters.npts(2));
        Gz=linspace(-0.5,0.5,parameters.npts(3));
        Gx=spdiags(Gx',0,parameters.npts(1),parameters.npts(1));
        Gy=spdiags(Gy',0,parameters.npts(2),parameters.npts(2));
        Gz=spdiags(Gz',0,parameters.npts(3),parameters.npts(3));
        Gx=parameters.dims(1)*kron(IdZ,kron(IdY,kron(Gx,Z)));
        Gy=parameters.dims(2)*kron(IdZ,kron(Gy,kron(IdX,Z)));
        Gz=parameters.dims(3)*kron(Gz,kron(IdY,kron(IdX,Z)));
        
    otherwise
        
        % Complain and bomb out
        error('incorrect number of spatial dimensions.');

end

% Project the infrastructure
H=kron(speye(prod(parameters.npts)),H);
K=kron(speye(prod(parameters.npts)),K);

% Preallocate relaxation matrix
R=spalloc(size(H,1),size(H,2),0);

% Build the relaxation matrix
if ~isempty(Ph)
    for n=1:numel(Ph)
        R=R+kron(spdiags(Ph{n}{2}(:),0,prod(parameters.npts),prod(parameters.npts)),Ph{n}{1});
    end
end
    
end

% You have my sympathies with rotations theory - this is one of the 
% most treacherous parts of spin dynamics. At the recent ENC, a guy
% approached Malcolm Levitt and declared that all rotation conventi-
% ons in his book are off by a sign. I could see Malcolm's face tur-
% ning white. Three hours later, when I completed my poster reading
% cycle and returned to Malcolm's poster, the two were still there,
% arguing heatedly over a laptop with Mathematica running. Somebody
% asked me what they were discussing. I said, "religion".
%
% IK's email to Fred Mentink-Vigier, 2014

