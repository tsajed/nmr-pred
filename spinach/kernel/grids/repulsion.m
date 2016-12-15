% Generates repulsion grids on a unit hypersphere. See the paper by
% Bak and Nielsen (http://dx.doi.org/10.1006/jmre.1996.1087) to get
% further information on the algorithm involved. Syntax:
%
%    [alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)
%
% Parameters:
%
%     npoints - number of points in the resulting spherical grid
%
%       ndims - hypersphere dimension: 2 returns a single-angle 
%               (beta) grid, 3 returns a two-angle grid (alpha,
%               beta), 4 returns a three-angle (alpha,beta,gam-
%               ma) spherical grid
%
%       niter - number of repulsion interations (simple clipped
%               gradient descent at the moment)
%
% Outputs:
%
%      alphas - alpha Euler angles of the grid, in radians,
%               zeros for single-angle grids
%
%       betas - beta Euler angles of the grid, in radians
%
%      gammas - gamma Euler angles of the grid, in radians,
%               zeros for two-angle grids
%
%     weights - point weights of the grid
% 
% Note: uniform weights are assigned at the moment, use the supp-
%       lied SHREWD function to generate optimal weights.
%
% i.kuprov@soton.ac.uk

function [alphas,betas,gammas,weights]=repulsion(npoints,ndims,niter)

% Check consistency
grumble(npoints,ndims,niter);

% Generate guess points
R=rand(ndims,npoints)-0.5;            

% Start the repulsion loop
for m=1:niter
    
    % Compute the forces
    F=zeros(size(R)); dots=R'*R;
    parfor n=1:npoints
        for k=1:npoints
            if n~=k
                distvec=(R(:,k)-R(:,n)); %#ok<PFBNS>
                distvec=distvec/norm(distvec);
                F(:,n)=F(:,n)+dots(n,k)*distvec;
            end
        end
    end
        
    % Apply the forces
    R_new=R-F/npoints;
    
    % Project points onto unit sphere
    R_new=R_new./kron(ones(ndims,1),sqrt(sum(R_new.^2,1))); 
    
    % Report the difference
    disp(['Iteration ' num2str(m) ', maximum point displacement: ' num2str(max(sqrt(sum((R-R_new).^2,2))))]);
    
    % Close the loop
    R=R_new;
    
end

% Get points
switch ndims
    
    case 2
        
        % In 2D case return polar angles
        [phi,~]=cart2pol(R(1,:),R(2,:));
        betas=phi'; alphas=0*betas; gammas=alphas;
        
    case 3
        
        % In 3D case return spherical angles
        [phi,theta,~]=cart2sph(R(1,:),R(2,:),R(3,:));
        betas=theta'+pi/2; alphas=phi'; gammas=0*alphas;
    
    case 4
        
        % In 4D case return Euler angles
        [gammas,betas,alphas]=quat2angle(R','ZYZ');
       
end

% Get weights
weights=ones(npoints,1)/npoints;

end

% Consistency enforcement
function grumble(npoints,ndims,niter)
if (~isnumeric(npoints))||(~isreal(npoints))||(numel(npoints)~=1)||(npoints<1)||(mod(npoints,1)~=0)
    error('npoints must be a positive integer.');
end
if (~isnumeric(niter))||(~isreal(niter))||(numel(niter)~=1)||(niter<1)||(mod(niter,1)~=0)
    error('niter must be a positive integer.');
end
if (~isnumeric(ndims))||(~isreal(ndims))||(numel(ndims)~=1)||(~ismember(ndims,[2 3 4]))
    error('ndims must be 2, 3 or 4.');
end
end
                                                          
% Let us beware of saying that there are laws in nature. There 
% are only necessities: there is nobody who commands, nobody
% who obeys, nobody who trespasses.
%
% Friedrich Nietzsche

