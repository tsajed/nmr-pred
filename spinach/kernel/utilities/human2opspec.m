% Converts user-friendly descriptions of spin states and operators into the 
% formal description (opspec) used by Spinach kernel. The function supports
% two types of calls:
%
% 1. If both inputs are strings, e.g.
%
%       [opspecs,coeffs]=human2opspec(spin_system,'Lz','13C')
%
% the function returns a list of single-spin opspecs for all spins with the
% specified name. In the example above, the list of Lz operator specificati-
% ons for all 13C nuclei in the system would be returned.
%
% Valid labels for operators in this type of call are 'E', 'Lz', 'L+', 'L-'
% and 'Tl,m'. In the latter case, l and m are integers. Valid labels for the
% spins are standard isotope names as well as 'electrons', 'nuclei', 'all'.
%
% 2. If the two inputs are a cell array of strings and a cell array of num-
% bers respectively, a product operator specification is produced, e.g.
% 
%       [opspecs,coeffs]=human2opspec(spin_system,{'Lz','L+'},{1,2})
%
% would return the LzL+ operator specification with Lz on spin 1 and L+ on
% spin 2. Valid labels for operators in the cell array are 'E', 'Lz', 'L+',
% 'L-' and 'Tl,m'. In the latter case, l and m are integers.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function [opspecs,coeffs]=human2opspec(spin_system,operators,spins)

% Check input validity
grumble(operators,spins);

% 'Lz','1H' type call returns a sum
if ischar(operators)&&ischar(spins)
    
    % Parse the specification
    switch spins
        
        case 'all'
            
            % Use all spins in the system
            spin_numbers=1:spin_system.comp.nspins;
            
        case 'electrons'
            
            % Count electrons
            spin_numbers=find(cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
            
        case 'nuclei'
            
            % Count non-electrons
            spin_numbers=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
            
        otherwise
            
            % Count specific spins
            spin_numbers=find(strcmp(spin_system.comp.isotopes,spins));
            
    end
    
    % Bomb out if the list of spins is empty
    if numel(spin_numbers)==0, error('no such spins in the system.'); end
    
    % Preallocate descriptors
    coeffs=zeros(numel(spin_numbers),1);
    opspecs=cell(numel(spin_numbers),1);
    
    % Use recursive calls
    for n=1:numel(spin_numbers)
        [opspecs(n),coeffs(n)]=human2opspec(spin_system,{operators},{spin_numbers(n)});
    end

    % Return the outcome
    return

% 'Lz',[1 2 4] type call returns a sum
elseif ischar(operators)&&isnumeric(spins)
    
    % Preallocate descriptors
    coeffs=zeros(numel(spins),1);
    opspecs=cell(numel(spins),1);
    
    % Use recursive calls
    for n=1:numel(spins)
        [opspecs(n),coeffs(n)]=human2opspec(spin_system,{operators},{spins(n)});
    end
    
    % Return the outcome
    return
    
% {'Lz','L+'},{1,4} type call returns an operator
elseif iscell(operators)&&iscell(spins)

    % Start with an empty opspec and a unit coefficient
    opspecs={zeros(1,spin_system.comp.nspins)}; coeffs=1;
    
    % Parse operator selection
    for n=1:numel(operators)
        
        switch operators{n}
            case 'E'
                
                % Unit operator
                opspecs{1}(spins{n})=0; coeffs=1*coeffs;
                
            case 'L+'
                
                % Raising operator
                opspecs{1}(spins{n})=1; coeffs=-sqrt(2)*coeffs;
                
            case 'Lz'
                
                % Z projection operator
                opspecs{1}(spins{n})=2; coeffs=1*coeffs;
                
            case 'L-'
                
                % Lowering operator
                opspecs{1}(spins{n})=3; coeffs=sqrt(2)*coeffs;
                
            otherwise
                
                % Validate the irreducible spherical tensor input
                if isempty(regexp(operators{n},'^T([\+\-]?\d+),([\+\-]?\d+)$','once'))
                    error('unrecognized operator specification.');
                end
                
                % Extract the quantum numbers
                indices=textscan(operators{n},'T%n,%n'); l=indices{1}; m=indices{2};
                
                % Validate the quantum numbers
                if (l<0)||(abs(m)>l)
                    error('invalid indices in irreducible spherical tensors.');
                end
                
                % Write the specification
                opspecs{1}(spins{n})=lm2lin(l,m); coeffs=1*coeffs;
                
        end
        
    end
    
else
    
    % Complain and bomb out
    error('invalid operator or state specification.');
    
end

end

% Input validation
function grumble(operators,spins)
if (~(ischar(operators)&&ischar(spins)))&&...
   (~(iscell(operators)&&iscell(spins)))&&...
   (~(ischar(operators)&&isnumeric(spins)))
    error('invalid operator or state specification.');
end
if iscell(operators)&&(~all(cellfun(@ischar,operators)))
    error('if operators is a cell array, all elements must be strings.');
end
if iscell(spins)&&(~all(cellfun(@isnumeric,spins)))
    error('if spins is a cell array, all elements must be numbers.');
end
end

% A "moral commandment" is a contradiction in terms. The moral is the
% chosen, not the forced; the understood, not the obeyed. The moral is
% the rational, and reason accepts no commandments.
%
% Ayn Rand

