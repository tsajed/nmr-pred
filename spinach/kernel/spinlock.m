% Analytical approximation to a spin locking process. This function oblite-
% rates all spin-spin correlations and all magnetization components other
% than those along the indicated direction. Parameters:
%
%      Lx         - X magnetization operator on the spins that
%                   should be locked
%
%      Ly         - Y magnetization operator on the spins that
%                   should be locked
%
%      rho        - state vector or a bookshelf stack thereof
%
%      direction  - direction in which the spins should be lo-
%                   cked, 'X' or 'Y'.
%
% The function destroys all magnetization except in the direction indicated,
% but does not prevent its subsequent evolution.
%
% i.kuprov@soton.ac.uk

function rho=spinlock(spin_system,Lx,Ly,rho,direction)

switch direction
    
    case 'X'
        
        % Destroy everything except for X magnetization
        rho=step(spin_system,Ly,rho,pi/2);
        rho=homospoil(spin_system,rho,'destroy');
        rho=step(spin_system,Ly,rho,-pi/2);
        
    case 'Y'

        % Destroy everything except for Y magnetization
        rho=step(spin_system,Lx,rho,pi/2);
        rho=homospoil(spin_system,rho,'destroy');
        rho=step(spin_system,Lx,rho,-pi/2);
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized spin locking direction.');
        
end

end

% To every dog - his bone and cage,
% To every wolf - his teeth and rage.
%
% Victor Tsoy

