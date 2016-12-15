% Plots grid integration quality as a function of spherical rank. The
% quality is defined as the norm of the residual of spherical harmon-
% ics or Wigner functions integrated using the grid provided. Syntax:
%
%  grid_profile=grid_test(alphas,betas,gammas,weights,max_rank,sfun)
%
% Parameters:
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
%    max_rank - maximum spherical rank to consider
%
%        sfun - spherical function type: for three-angle
%               grids use 'D_lmn', for two-angle grids use
%               'Y_lm', for single-angle grids use 'Y_l0'.
%
% The output is a vector of residual norms in each spherical rank.
%
% i.kuprov@soton.ac.uk

function grid_profile=grid_test(alphas,betas,gammas,weights,max_rank,sfun)

% Preallocate the answer
grid_profile=zeros(max_rank+1,1);

% Loop over spherical ranks
for l=0:max_rank
    
    % Preallocate Wigner matrix
    D=zeros(2*l+1);
    
    % Loop over grid points
    parfor n=1:numel(alphas)
        D=D+weights(n)*wigner(l,alphas(n),betas(n),gammas(n));
    end
    
    % Update grid profile
    if strcmp(sfun,'D_lmn')
        grid_profile(l+1)=norm(D)-krondelta(0,l);
    elseif strcmp(sfun,'Y_lm')
        grid_profile(l+1)=norm(D(l+1,:))-krondelta(0,l);
    elseif strcmp(sfun,'Y_l0')
        grid_profile(l+1)=norm(D(l+1,l+1))-krondelta(0,l);
    else
        error('unknown diagnostics function.');
    end
    
    % Update the user
    disp(['Spherical rank ' num2str(l) ', residual ' sfun ' norm: ' num2str(grid_profile(l+1))]);
    
end
    
% Do the plotting
plot(0:max_rank,abs(grid_profile));
xlabel('spherical rank');
ylabel('integration residual');

end

% Art is what you can get away with.
%
% Andy Warhol

