% Optimal step for time propagation under a given Hamiltonian or a given
% Liouvillian superoperator. The function uses the 1-norm (which is the
% cheapest variety and the safest one in terms of matrix scaling).
%
% Syntax:          [timestep,nsteps]=stepsize(L,interval)
%
% Parameters:
%
%     L   -      Liouvillian operator or superoperator inder which the 
%                evolution is taking place.
%
%     interval - the total evolution time, must be specified if the op-
%                timal number of steps is also wanted in the output.
%
% Outputs:
%
%     timestep - optimal time step
%
%     nsteps   - number of steps covering the interval, if the interval
%                was specified
%
% Note: if a zero Liouvillian is supplied, the function would return the
%       following values: timestep=0, nsteps=1
%
% i.kuprov@soton.ac.uk

function [timestep,nsteps]=stepsize(L,interval)

% Catch incorrect calls
if (nargin==1)&&(nargout==2)
    error('time interval is required to compute the number of steps.');
elseif (~isnumeric(L))
    error('L must be numeric.');
elseif (nargin==2)&&((~isnumeric(interval))||(~isreal(interval))||(~isscalar(interval)))
    error('interval muust be a real number.');
end

% Estimate the scaling coefficient   
scaling=max([1 norm(L,1)]);

% Adapt to the output style
switch nargout
    
    case 1
        
        % Get the optimal time step
        timestep=1/scaling;
    
    case 2
               
        % Get the optimal integer number of steps
        nsteps=ceil(abs(interval)*scaling);
        
        % Get the time step
        timestep=interval/nsteps;
        
end

end

% Whenever anyone accuses some person of being 'unfeeling' he means that
% that person is just. He means that that person has no causeless emotions
% and will not grant him a feeling which he does not deserve... justice is
% the opposite of charity.
%
% Ayn Rand, "Atlas Shrugged"

