% Computes Gauss-Legendre points and weights in [a,b] interval
% with accuracy order n.
%
% i.kuprov@soton.ac.uk

function [x,w]=gaussleg(a,b,n)

% Initial guess for the nodes in [-1 1]
x=cos((2*(0:n)'+1)*pi/(2*n+2))+(0.27/(n+1))*sin(pi*linspace(-1,1,n+1)'*n/(n+2));

% Newton-Raphson refinement
V=zeros(n+1,n+2); prev_x=0;
while max(abs(x-prev_x))>eps
    V(:,1)=1; V(:,2)=x;
    for k=2:(n+1)
        V(:,k+1)=((2*k-1)*x.*V(:,k)-(k-1)*V(:,k-1))/k;
    end
    dV=(n+2)*(V(:,n+1)-x.*V(:,n+2))./(1-x.^2);
    prev_x=x; x=x-V(:,n+2)./dV;
end

% Compute the weights
w=(b-a)./((1-x.^2).*dV.^2)*((n+2)/(n+1))^2;

% Map from [-1,1] to [a,b]
x=(a*(1-x)+b*(1+x))/2;

end

% Life is a fountain of delight; but where the rabble also 
% drinks all wells are poisoned.
%
% Friedrich Nietzsche

