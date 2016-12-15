% Simpson's numerical integration.
%
% Z = SIMPS(Y) computes an approximation of the integral of Y using 
% Simpson's method (with unit spacing). To compute the integral for
% spacing different from one, multiply Z by the spacing increment.
%
% For vectors, SIMPS(Y) is the integral of Y. For matrices, SIMPS(Y)
% is a row vector with the integral over each column. For N-D arrays,
% SIMPS(Y) works across the first non-singleton dimension.
%
% Z = SIMPS(X,Y) computes the integral of Y with respect to X using 
% Simpson's rule.  X and Y must be vectors of the same length, or X
% must be a column vector and Y an array whose first non-singleton 
% dimension is length(X). SIMPS operates along this dimension.
%
% Z = SIMPS(X,Y,DIM) or SIMPS(Y,DIM) integrates across dimension DIM
% of Y. The length of X must be the same as size(Y,DIM).
%
% damien.garcia@biomecardio.com
% i.kuprov@soton.ac.uk

function z = simps(x,y,dim)

% Make sure x and y are column vectors, or y is a matrix.
perm = []; nshifts = 0;
if nargin == 3 % simps(x,y,dim)
  perm = [dim:max(ndims(y),dim) 1:dim-1];
  yp = permute(y,perm);
  [m,n] = size(yp);
elseif nargin==2 && isscalar(y) % simps(y,dim)
  dim = y; y = x;
  perm = [dim:max(ndims(y),dim) 1:dim-1];
  yp = permute(y,perm);
  [m,n] = size(yp);
  x = 1:m;
else % simps(y) or simps(x,y)
  if nargin < 2, y = x; end
  [yp,nshifts] = shiftdim(y);
  [m,n] = size(yp);
  if nargin < 2, x = 1:m; end
end
x = x(:);
if length(x) ~= m
  if isempty(perm) % dim argument not given
    error('MATLAB:simps:LengthXmismatchY',...
          'LENGTH(X) must equal the length of the first non-singleton dimension of Y.');
  else
    error('MATLAB:simps:LengthXmismatchY',...
          'LENGTH(X) must equal the length of the DIM''th dimension of Y.');
  end
end

% The output size for [] is a special case when DIM is not given.
if isempty(perm) && isequal(y,[])
  z = zeros(1,class(y));
  return;
end

% Use TRAPZ if m<3
if m<3
    if exist('dim','var')
        z = trapz(x,y,dim);
    else
        z = trapz(x,y);
    end
    return
end

% Simpson's rule
y = yp;
clear yp

dx = repmat(diff(x,1,1),1,n);
dx1 = dx(1:end-1,:);
dx2 = dx(2:end,:);

alpha = (dx1+dx2)./dx1/6;
a0 = alpha.*(2*dx1-dx2);
a1 = alpha.*(dx1+dx2).^2./dx2;
a2 = alpha.*dx1./dx2.*(2*dx2-dx1);

z = sum(a0(1:2:end,:).*y(1:2:m-2,:) +...
    a1(1:2:end,:).*y(2:2:m-1,:) +...
    a2(1:2:end,:).*y(3:2:m,:),1);

if rem(m,2) == 0 % Adjusting if length(x) is even
    state0 = warning('query','MATLAB:nearlySingularMatrix');
    state0 = state0.state;
    warning('off','MATLAB:nearlySingularMatrix')
    C = vander(x(end-2:end))\y(end-2:end,:);
    z = z + C(1,:).*(x(end,:).^3-x(end-1,:).^3)/3 +...
        C(2,:).*(x(end,:).^2-x(end-1,:).^2)/2 +...
        C(3,:).*dx(end,:);
    warning(state0,'MATLAB:nearlySingularMatrix')
end

% Resizing
siz = size(y); siz(1) = 1;
z = reshape(z,[ones(1,nshifts),siz]);
if ~isempty(perm), z = ipermute(z,perm); end

end


