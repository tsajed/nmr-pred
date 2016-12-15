% Sampling function for ZFS parameter distributions.
%
% <http://spindynamics.org/wiki/index.php?title=Zfs_sampling.m>

function [D,E,W]=zfs_sampling(npoints_d,npoints_e,tol)

% Check consistency
grumble(npoints_d,npoints_e,tol);

% Generate Gauss-Legendre point set for D/D1
[X,WX]=lgwt(npoints_d,-2,2);

% Get the standard deviation for unit FWHM
sigma=1/(2*sqrt(2*log(2)));

% Refract weights through a double Gaussian
WX=WX.*(normpdf(X,-1,sigma)+normpdf(X,+1,sigma)); WX=WX/sum(WX);

% Plot the double Gaussian
subplot(1,2,1); plot(X,normpdf(X,-1,sigma)+normpdf(X,+1,sigma),'r-');
title('D/D0 distribution'); xlabel('D/D0');
ylabel('probability density'); axis tight;

% Generate Gauss-legendre point set for E/D
[Y,WY]=lgwt(npoints_e,0,1/3);

% Refract weights through a quadratic function
WY=WY.*(-(Y-0.25).^2+0.0625); WY=WY/sum(WY);

% Plot the quadratic function
subplot(1,2,2); plot(Y,-(Y-0.25).^2+0.0625,'r-');
title('E/D distribution'); xlabel('E/D');
ylabel('probability density'); axis tight;

% Kron the weights
D=kron(X,ones(size(Y)));
E=D.*kron(ones(size(X)),Y);
W=kron(WX,WY); W=W/sum(W);

% Ignore small weights
D(W<tol)=[]; E(W<tol)=[]; W(W<tol)=[];

end

% Consistency enforcement
function grumble(npoints_d,npoints_e,tol)
if (~isnumeric(npoints_d))||(~isreal(npoints_d))||...
   (npoints_d<5)||(mod(npoints_d,1)~=0)
    error('npoints_d must be a real integer greater than 5.');
end
if (~isnumeric(npoints_e))||(~isreal(npoints_e))||...
   (npoints_e<5)||(mod(npoints_e,1)~=0)
    error('npoints_e must be a real integer greater than 5.');
end
if (~isnumeric(tol))||(~isreal(tol))||(tol<0)
    error('tol must be a positive real number much smaller than 1.');
end
end

% It is a cliche that most cliches are true, but 
% then, like most cliches, that cliche is untrue.
%
% Stephen Fry

