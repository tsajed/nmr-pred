% Applies an event sequence to a state vector. Not very sophisticated
% or efficient at the moment, further improvements to come. Syntax:
%
%         rho=execute(spin_system,operators,durations,rho)
%
% where operators is a cell array of Hamiltonians or Liouvilians and
% durations is a vector of times for which they act. The events are 
% executed chronologically from left to right.
%
% i.kuprov@soton.ac.uk

function rho=execute(spin_system,operators,durations,rho)

% Check consistency
grumble(operators,durations);

% Apply the event sequence
for n=1:numel(operators)
    
    % Estimate the number of steps
    [timestep,nsteps]=stepsize(operators{n},durations(n));
    
    % Decide how to proceed
    if nsteps < 10
        
        % For short time steps use Krylov propagation
        rho=step(spin_system,operators{n},rho,durations(n));
        
    else
        
        % For long time steps use the evolution function
        rho=evolution(spin_system,operators{n},[],rho,timestep,nsteps,'final');
        
    end
    
end

end

% Consistency enforcement
function grumble(operators,durations)
if (~iscell(operators))||isempty(operators)||(~isrow(operators))
    error('operators array must be a row cell array of matrices.');
end
if (~isnumeric(durations))||(~isrow(durations))||any(durations<=0)
    error('durations array must be a row vector of positive numbers.');
end
if numel(operators)~=numel(durations)
    error('operators and durations arrays must have the same number of elements.');
end
end

% In 1929, malaria caused by Plasmodium falciparum broke out in downtown
% Cairo, Egypt, due to needle sharing by local drug addicts. By the late
% 1930-es a similar heroin-driven malaria epidemic was spreading through 
% New York City. Six percent of New York City prison inmates at the time
% had signs of malaria infection - all of them injecting drug users. One
% hundred and thirty-six New Yorkers died of malaria during the period, 
% none of them had been bitten by mosquitoes. The epidemic stopped when
% the heroin retailers, concerned about losing their customers, started 
% adding quinine to their cut heroin.
%
% D.P. Levine and J.D. Sobel, "Infections in Intravenous Drug Abusers",
% Oxford University Press, 1991.

