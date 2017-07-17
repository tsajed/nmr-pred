% Returns a finite-difference representation of the Laplacian for a 3D
% array with a user-specified finite difference stencil size. The re-
% sulting operator is a sparse matrix designed to act on the vectori-
% zation of the 3D array. The dimensions of the 3D array are assumed
% to be ordered as [X Y Z]. Syntax:
%
%                    L=fdlap(npoints,extents,nstenc)
%
% The following parameters are needed:
% 
%     dims    -  a three-element vector specifying the number of
%                discretization points in each dimension of the
%                3D cube of data that the operator will be acting
%                on, ordered as [X Y Z].
%
%     extents -  a three-element vector specifying axis extents,
%                ordered as [X Y Z].
%
%     nstenc  -  number of finite-difference stencil points for
%                the finite-difference approximations.
%
% The resulting operator is a sparse matrix designed to act on the vec- 
% torization of rho. The dimensions of rho are assumed to be ordered 
% as [X Y Z].
%
% Periodic boundary condition is used.
%
% i.kuprov@soton.ac.uk

function L=fdlap(dims,extents,nstenc)

% Check consistency
grumble(dims,extents,nstenc);

switch numel(dims)
    
    case 1
        
        % Get differentiation matrices
        Dxx=fdmat(dims(1),nstenc,2);
        
        % Normalize differentiation matrices
        Dxx=(dims(1)/extents(1))^2*Dxx;
        
        % Compute the Laplacian
        L=Dxx;
    
    case 2
        
        % Get differentiation matrices
        Dxx=fdmat(dims(1),nstenc,2);
        Dyy=fdmat(dims(2),nstenc,2);
        
        % Normalize differentiation matrices
        Dxx=(dims(1)/extents(1))^2*Dxx;
        Dyy=(dims(2)/extents(2))^2*Dyy;
        
        % Compute the Laplacian
        L=kron(Dyy,speye(dims(1)))+...
          kron(speye(dims(2)),Dxx);
    
    case 3

        % Get differentiation matrices
        Dxx=fdmat(dims(1),nstenc,2);
        Dyy=fdmat(dims(2),nstenc,2);
        Dzz=fdmat(dims(3),nstenc,2);
        
        % Normalize differentiation matrices
        Dxx=(dims(1)/extents(1))^2*Dxx;
        Dyy=(dims(2)/extents(2))^2*Dyy;
        Dzz=(dims(3)/extents(3))^2*Dzz;
        
        % Compute the Laplacian
        L=kron(kron(Dzz,speye(dims(2))),speye(dims(1)))+...
          kron(kron(speye(dims(3)),Dyy),speye(dims(1)))+...
          kron(kron(speye(dims(3)),speye(dims(2))),Dxx);

    otherwise
        
        % Complain and bomb out
        error('incorrect number of spatial dimensions.');
        
end

end

% Consistency enforcement
function grumble(dims,extents,nstenc)
if (~isnumeric(dims))||(~isreal(dims))||(any(dims<1))||any(mod(dims,1)~=0)
    error('npoints must be a three-element vector of positive integers.');
end
if (~isnumeric(extents))||(~isreal(extents))||(any(extents<=0))
    error('extents must be an array of positive real numbers.');
end
if any(dims<nstenc)
    error('array dimension is not big enough for the finite difference stencil specified.');
end
if (mod(nstenc,1)~=0)||(mod(nstenc,2)~=1)||(nstenc<3)
    error('the number of stencil points must be an odd integer greater or equal to 3.');
end
end

% We owe no morality to those who hold us under a gun.
% 
% Ayn Rand
