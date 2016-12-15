% Computes SHREWD weights for a given two- or three-angle spherical
% grid. See the paper by Eden and Levitt for details on now the al-
% gorithm works: http://dx.doi.org/10.1006/jmre.1998.1427 Syntax:
%
%      weights=shrewd(alphas,betas,gammas,max_rank,max_error)
%
% Parameters:
%
%      alphas - alpha Euler angles of the grid, in radians
%
%       betas - beta Euler angles of the grid, in radians
%
%      gammas - gamma Euler angles of the grid,in radians,
%               set to all-zeros for two-angle grids
%
%    max_rank - maximum spherical rank to take into consi-
%               deration when minimizing residuals
%
%   max_error - maximum residual absolute error per spheri-
%               cal function
%
% The output is a vector of grid weights for each [alpha beta gamma]
% point supplied.
%
% Note: for a given arrangement of angles, this is the most consis-
%       tent weight selection procedure in the literature.
%
% i.kuprov@soton.ac.uk

function weights=shrewd(alphas,betas,gammas,max_rank,max_error)

% Decide the grid type
if all(gammas==0)

    % Preallocate spherical harmonic matrix
    H=zeros(lm2lin(max_rank,-max_rank)+1,numel(alphas));

    % Fill spherical harmonic matrix
    for k=1:numel(alphas)
        for l=0:max_rank
            D=wigner(l,alphas(k),betas(k),gammas(k));
            for m=l:-1:-l
                H(lm2lin(l,m)+1,k)=D(l+1,l+m+1);
            end
        end
    end

    % Get the right hand side vector
    v=max_error*ones(lm2lin(max_rank,-max_rank)+1,1);
    v(1)=1-max_error;
    
else
    
    % Preallocate Wigner function matrix
    H=zeros(lmn2lin(max_rank,-max_rank,-max_rank),numel(alphas));

    % Fill Wigner function matrix
    for k=1:numel(alphas)
        for l=0:max_rank
            D=wigner(l,alphas(k),betas(k),gammas(k));
            for m=l:-1:-l
                for n=l:-1:-l
                    H(lmn2lin(l,m,n),k)=D(l+m+1,l+n+1);
                end
            end
        end
    end

    % Get the right hand side vector
    v=max_error*ones(lmn2lin(max_rank,-max_rank,-max_rank),1);
    v(1)=1-max_error;
    
end

% Compute the weights
weights=real(H\v);
weights=weights/sum(weights);

end

% A good theory explains a lot but postulates little.
%
% Richard Dawkins

