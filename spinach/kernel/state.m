% State operators (Hilbert space) and state vectors (Liouville space).
%
% <http://spindynamics.org/wiki/index.php?title=State.m>

function rho=state(spin_system,states,spins,method)

% Validate the input
grumble(spin_system,states,spins);

% Default is to use consistent state norms
if ~exist('method','var'), method='exact'; end

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Choose the state vector generation methos
        switch method
            
            case 'cheap'
                
                % Parse the specification
                [opspecs,coeffs]=human2opspec(spin_system,states,spins);
                
                % Compute correlation orders
                correlation_orders=sum(logical(spin_system.bas.basis),2);
                
                % Get the state vector by searching the state list
                rho=spalloc(size(spin_system.bas.basis,1),1,numel(opspecs));
                parfor n=1:numel(opspecs) %#ok<*PFBNS>
                    
                    % Find states with the same correlation order
                    possibilities=(correlation_orders==nnz(opspecs{n}));
                    
                    % Pin down the required state
                    for k=find(opspecs{n})
                        possibilities=and(possibilities,spin_system.bas.basis(:,k)==opspecs{n}(k)); 
                    end

                    % Double-check
                    if nnz(possibilities)>1
                        error('basis descriptor ambiguity detected.');
                    elseif nnz(possibilities)<1
                        error('the requested state is not present in the basis.');
                    end
                    
                    % Add to the total
                    rho=rho+coeffs(n)*sparse(possibilities);
                    
                end
                
            case 'exact'
                
                % Apply a left side product superoperator to the unit state
                rho=operator(spin_system,states,spins,'left')*unit_state(spin_system);
                
            case 'chem'
                
                % Parse the specification
                [opspecs,coeffs]=human2opspec(spin_system,states,spins);
                
                % Preallocate the state vector
                rho=spalloc(size(spin_system.bas.basis,1),1,0);
                
                % Sum the states with concentrations
                for n=1:numel(opspecs)
                    
                    % Identify active spins
                    active_spins=find(opspecs{n});
                    
                    % Find out which chemical species they are in
                    species=true(1,numel(spin_system.chem.parts)); 
                    for k=1:numel(active_spins)
                        species=species&cellfun(@(x)ismember(active_spins(k),x),spin_system.chem.parts);
                    end
                    
                    % Check state validity
                    if nnz(species)~=1
                        error('the spin state requested crosses chemical species boundaries.');
                    end
                    
                    % Adjust the coefficient
                    coeffs(n)=coeffs(n)*spin_system.chem.concs(species);
                    
                    % Get the state vector
                    rho=rho+coeffs(n)*p_superop(spin_system,opspecs{n},'left')*unit_state(spin_system);
                    
                end
                
            otherwise
                
                % Complain and bomb out
                error('unknown state generation method.');
                
        end
        
    case 'zeeman-liouv'

        % Apply a left side product superoperator to the unit state
        rho=operator(spin_system,states,spins,'left')*unit_state(spin_system); 
                
    case 'zeeman-hilb'
        
        % Generate a Hilbert space operator
        rho=operator(spin_system,states,spins);
        
        % Match normalization to Liouville space
        rho=rho./sqrt(size(rho,1));
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
    
end

end

% Input validation function
function grumble(spin_system,states,spins)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~(ischar(states)&&ischar(spins)))&&...
   (~(iscell(states)&&iscell(spins)))&&...
   (~(ischar(states)&&isnumeric(spins)))
    error('invalid state specification.');
end
if iscell(states)&&iscell(spins)&&(numel(states)~=numel(spins))
    error('spins and operators cell arrays should have the same number of elements.');
end
if iscell(states)&&any(~cellfun(@ischar,states))
    error('all elements of the operators cell array should be strings.');
end
if iscell(spins)&&any(~cellfun(@isnumeric,spins))
    error('all elements of the spins cell array should be integer numbers.');
end
end

% It worked.
%
% J. Robert Oppenheimer (after witnessing the first atomic detonation) 

