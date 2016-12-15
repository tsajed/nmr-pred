% Returns arbitrary-order central finite-difference differentiation 
% matrix (sparse) with unit spacing and periodic boundary conditions.
% Syntax:
%
%                      D=fdmat(dim,npoints,order)
%
% Parameters:
%
%     dim       - dimension of the column vector to be
%                 differentiated
%
%     nstenc    - number of points in the finite diffe-
%                 rence stencil
%
%     order     - order of the derivative required
%
% i.kuprov@soton.ac.uk

function D=fdmat(dim,nstenc,order)

% Check consistency
grumble(dim,nstenc,order);

% Preallocate the answer
D=spalloc(dim,dim,dim*nstenc);

% Wraparound fill with centered schemes
stencil=((1-nstenc)/2):((nstenc-1)/2);
w=fdweights(0,stencil,order);
for n=1:dim
    D(n,mod(stencil+n-1,dim)+1)=w(end,:); %#ok<SPRIX>
end

end

% Consistency enforcement
function grumble(dim,nstenc,order)
if (dim<1)||(nstenc<1)||(order<1)||(mod(dim,1)~=0)||(mod(nstenc,1)~=0)||(mod(order,1)~=0)
    error('all input parameters must be positive integers.');
end
if dim<3
    error('minimum differentiation matrix dimension is 3.');
end
if mod(nstenc,2)~=1
    error('the number of stencil points must be odd.');
end
if order>=nstenc
    error('derivative order must be smaller than the stencil size.');
end
end

% We may eventually come to realize that chastity 
% is no more a virtue than malnutrition.
%
% Alex Comfort

