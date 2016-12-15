% Trajectory analysis function. Plots the time dependence of the densi-
% ty matrix norm, partitioned into user-specified property classes.
%
% Call syntax:
%
%     trajan(spin_system,trajectory,property)
%
% Arguments:
%
%     trajectory - a stack of state vectors of any length. The
%                  number of rows in the trajectory array must
%                  match the number of states in the basis.
%
%     property   - If set to 'correlation_order', returns the
%                  time dependence of the total populations of
%                  one-spin, two-spin, three-spin, etc. corre-
%                  lations in the system.
%
%                  If set to 'coherence_order', returns the ti-
%                  me dependence of different orders of coheren-
%                  ce in the system, where a coherence order is
%                  defined as the sum of projection quantum num-
%                  bers in the spherical tensor representation
%                  of each state.
%
%                  If set to 'total_each_spin', returns the ti-
%                  me dependence of total state space populati-
%                  on that involves each individual spin in the
%                  system in any way (all local populations and
%                  coherences of the spin as well as all of its 
%                  correlations to all third party spins).
%
%                  If set to 'local_each_spin', returns the ti-
%                  me dependence of the population of the sub-
%                  space of states that are local to each indi-
%                  vidual spin and do not involve any correla-
%                  tions to other spins in the system.
%
%                  If set to 'level_populations', returns the
%                  populations of the Zeeman energy levels. 
%
% The trajectory would usually come out of the evolution() run from a
% given starting point under a given Liouvillian.
%
% Note: this function is only applicable to the trajectories recorded
%       in spherical tensor basis sets.
%
% gareth.charnock@oerc.ox.ac.uk
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function trajan(spin_system,trajectory,property)

% Check the input
grumble(spin_system,trajectory);

% Set the defaults
if ~exist('property','var'), property='correlation_order'; end

% Determine how to proceed
switch property
    
    case 'correlation_order'
        
        % Determine the correlation order of each state
        correlation_orders=sum(logical(spin_system.bas.basis),2);
        
        % Find out which correlation orders are present
        unique_correlation_orders=unique(correlation_orders);
        
        % Preallocate the norm trajectory array
        result=zeros(numel(unique_correlation_orders),size(trajectory,2));
        
        % Loop over the unique correlation orders that are present
        for n=1:numel(unique_correlation_orders)
            
            % Find the subspace of the given correlation order
            subspace_mask=(correlation_orders==unique_correlation_orders(n));
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=trajectory(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create a legend
        legend_text=cell(numel(unique_correlation_orders),1);
        for n=1:numel(unique_correlation_orders)
            legend_text{n}=[num2str(unique_correlation_orders(n)) '-spin order'];
        end
        
        % Create an axis label
        label_text='correlation order amplitude';
        
    case 'coherence_order'
        
        % Determine projection quantum numbers of the basis
        [~,M]=lin2lm(spin_system.bas.basis);
        
        % Determine the coherence order of each state
        coherence_orders=sum(M,2);
        
        % Find out which coherence orders are present
        unique_coherence_orders=unique(coherence_orders);
        
        % Preallocate the norm trajectory array
        result=zeros(numel(unique_coherence_orders),size(trajectory,2));
        
        % Loop over the unique coherence orders that are present
        for n=1:numel(unique_coherence_orders)
            
            % Find the subspace of the given coherence order
            subspace_mask=(coherence_orders==unique_coherence_orders(n));
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=trajectory(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create a legend
        legend_text=cell(numel(unique_coherence_orders),1);
        for n=1:numel(unique_coherence_orders)
            legend_text{n}=['coherence order ' num2str(unique_coherence_orders(n))];
        end
        
        % Create an axis label
        label_text='coherence order amplitude';
        
    case 'total_each_spin'
        
        % Preallocate the norm trajectory array
        result=zeros(spin_system.comp.nspins,size(trajectory,2));
        
        % Loop over spins in the system
        for n=1:spin_system.comp.nspins
            
            % Find the subspace of states that involve the current spin
            subspace_mask=(spin_system.bas.basis(:,n)~=0);
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=trajectory(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create an axis label
        label_text='density touching each spin';
        
    case 'local_each_spin'
        
        % Preallocate the norm trajectory array
        result=zeros(spin_system.comp.nspins,size(trajectory,2));
        
        % Loop over spins in the system
        for n=1:spin_system.comp.nspins
            
            % Find the subspace of states that are local to current spin
            subspace_mask=(spin_system.bas.basis(:,n)~=0)&...
                          (sum(spin_system.bas.basis,2)==spin_system.bas.basis(:,n));
            
            % Get the part of the trajectory belonging to the subspace
            subspace_trajectory=trajectory(subspace_mask,:);
            
            % Get the norm of the trajectory
            result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1));
            
        end
        
        % Create an axis label
        label_text='density local to each spin';
        
    case 'level_populations'
        
        % Move trajectory into the Zeeman basis set
        trajectory=sphten2zeeman(spin_system)*trajectory;
        
        % Find out the number of energy levels
        nlevels=sqrt(size(trajectory,1));
        
        % Preallocate population dynamics array
        result=zeros(nlevels,size(trajectory,2));
        
        % Extract the populations
        for n=1:size(trajectory,2)
            result(:,n)=real(diag(reshape(trajectory(:,n),[nlevels nlevels])));
        end
        
        % Create an axis label
        label_text='energy level populations';
        
    otherwise
        
        % Complain and bomb out
        error('unknown property.');
        
end

% Do the plotting
plot(result'); xlabel('time / points'); 
ylabel(label_text); axis tight;
set(gca,'yscale','log'); 
if exist('legend_text','var')
    legend(legend_text,'Location','NorthEast');
end

end

% Consistency enforcement
function grumble(spin_system,trajectory)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('trajectory analysis is only available for sphten-liouv formalism.');
end
if ~isnumeric(trajectory)
    error('trajectory should be an array of doubles.');
end
if size(trajectory,1)~=size(spin_system.bas.basis,1)
    error('trajectory dimension should match basis dimension.');
end
end

% And this is the whole shabby secret: to some men, the sight of
% an achievement is a reproach, a reminder that their own lives 
% are irrational, and that there is no loophole - no escape from
% reason and reality.
%
% Ayn Rand

