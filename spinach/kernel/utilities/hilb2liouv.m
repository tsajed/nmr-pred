% Converts Hilbert space operators into Liouville space superoperators
% or state vectors. Syntax:
%
%                           L=hilb2liouv(H,conv_type)
%
% Parameters:
%
%          H - a Hilbert space operator
%
%  conv_type - the type of Liouville space superoperator to return:
%
%                  'left' - produces left side product superoperator
%
%                 'right' - produces right side product superoperator
%
%                  'comm' - produces commutation superoperator
%
%                 'acomm' - produces anticommutation superoperator
%
%              'statevec' - stretches the operator into a state vector
%
% i.kuprov@soton.ac.uk

function L=hilb2liouv(H,conv_type)

% Check consistency
grumble(H,conv_type);

% Prepare a unit matrix
unit=speye(size(H));

% Decide how to proceed
switch conv_type
    
    case 'comm'
        
        % Compute a commutation superoperator
        L=kron(unit,H)-kron(transpose(H),unit);
  
    case 'acomm'
        
        % Compute an anticommutation superoperator
        L=kron(unit,H)+kron(transpose(H),unit);
        
    case 'left'
        
        % Compute a left side product superoperator
        L=kron(unit,H);
        
    case 'right'
        
        % Compute a right side product superoperator
        L=kron(transpose(H),unit);
        
    case 'statevec'
        
        % Stretch into a state vector
        L=H(:);
        
    otherwise
        
        % Complain and bomb out
        error('unknown conversion type specification.');
        
end

end

% Consistency enforcement
function grumble(H,conv_type)
if ~isnumeric(H)
    error('H parameter must be numeric.');
end
if ~ischar(conv_type)
    error('conv_type parameter must be a character string.');
end
end

% If it wasn't for the fun and money, I really don't know 
% why I'd bother.
%
% Terry Pratchett

