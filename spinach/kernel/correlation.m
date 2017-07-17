% Correlation order selection function -- keeps only the specified orders
% of spin correlation in the state vector. Syntax:
%
%       rho=correlation(spin_system,rho,correlation_orders,spins)
%
% Arguments:
%
%   rho                -  a state vector or a horizontal stack thereof
%
%   correlation_orders -  a row vector of correlation orders to keep
%
%   spins              -  which spins to count (e.g. '1H', '13C', 'all')
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function rho=correlation(spin_system,rho,correlation_orders,spins)

% Set the default to all spins
if ~exist('spins','var'), spins='all'; end

% Check consistency
grumble(spin_system,rho,correlation_orders)

% Store dimension statistics
spn_dim=size(spin_system.bas.basis,1);
spc_dim=numel(rho)/spn_dim;
problem_dims=size(rho);

% Fold indirect dimensions
rho=reshape(rho,[spn_dim spc_dim]);

% Parse the spin specification
if ~isnumeric(spins)
    if strcmp(spins,'all')
        spins=1:length(spin_system.comp.isotopes);
    else
        spins=find(strcmp(spins,spin_system.comp.isotopes));
    end
end

% Compute the order of correlation for each basis state
correlation_orders_present=sum(logical(spin_system.bas.basis(:,spins)),2);

% Wipe all correlation orders except those specified by the user
state_mask=zeros(size(spin_system.bas.basis,1),1);
for n=correlation_orders
    state_mask=state_mask|(correlation_orders_present==n);
end
    
% Apply the mask
rho(~state_mask,:)=0;

% Unfold indirect dimensions
rho=reshape(rho,problem_dims);

% Report overly destructive calls
if norm(rho,1)<1e-10
    report(spin_system,'WARNING - all magnetization appears to have been destroyed by this call.');
end
    
end

% Consistency enforcement
function grumble(spin_system,rho,correlation_orders)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('analytical correlation order selection is only available for sphten-liouv formalism.');
end
if ~isnumeric(rho)
    error('the state vector(s) must be numeric.');
end
if (~isnumeric(correlation_orders))||(~isvector(correlation_orders))||...
     any(correlation_orders<0)||any(mod(correlation_orders,1)~=0)
    error('correlation_orders parameter must be a vector of non-negative integers.');
end
end

% The first iteration of the Spin Dynamics course (http://spindynamics.org) 
% was so difficult that every single student has dropped out by about Lecture
% 10. The rest of the course was read to Rusty, Dusty, Scratchy, Patchy and
% Scruffy, the five plastic chairs in IK's Durham office - they made for an 
% excellent (if a bit shy so far as questions were concerned) audience. To 
% this day, the total number of students who verifiably understood the whole
% of that course can be counted on the fingers of one hand.
