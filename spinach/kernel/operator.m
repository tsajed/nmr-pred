% Hilbert space operators and Liouville space superoperators.
%
% <http://spindynamics.org/wiki/index.php?title=Operator.m>

function A=operator(spin_system,operators,spins,operator_type)

% Validate the input
grumble(spin_system,operators,spins);

% The default type is commutation superoperator
if ~exist('operator_type','var'), operator_type='comm'; end

% Parse the specification
[opspecs,coeffs]=human2opspec(spin_system,operators,spins);

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Get the Liouville space superoperator
        local_answers=cell(numel(opspecs),1);
        parfor n=1:numel(opspecs)
            local_answers{n}=coeffs(n)*p_superop(spin_system,opspecs{n},operator_type);
        end
        A=spalloc(size(local_answers{1}(1),1),...
                  size(local_answers{1}(2),2),...
                  sum(cellfun(@nnz,local_answers)));
        for n=1:numel(opspecs)
            A=A+local_answers{n};
        end
                
    case 'zeeman-hilb'
        
        % Get the Hilbert space operator
        local_answers=cell(numel(opspecs),1);
        parfor n=1:numel(opspecs)
            B=sparse(coeffs(n));
            for k=1:numel(opspecs{n})
                [L,M]=lin2lm(opspecs{n}(k));
                ist=irr_sph_ten(spin_system.comp.mults(k),L); %#ok<PFBNS>
                B=kron(B,ist{L-M+1});
            end
            local_answers{n}=B;
        end
        A=spalloc(size(local_answers{1}(1),1),...
                  size(local_answers{1}(2),2),...
                  sum(cellfun(@nnz,local_answers)));
        for n=1:numel(opspecs)
            A=A+local_answers{n};
        end

    case 'zeeman-liouv'
        
        % Get the Liouville space superoperator
        local_answers=cell(numel(opspecs),1);
        parfor n=1:numel(opspecs)
            B=sparse(coeffs(n));
            for k=1:numel(opspecs{n})
                [L,M]=lin2lm(opspecs{n}(k));
                ist=irr_sph_ten(spin_system.comp.mults(k),L); %#ok<PFBNS>
                B=kron(B,ist{L-M+1});
            end
            local_answers{n}=hilb2liouv(B,operator_type);
        end
        A=spalloc(size(local_answers{1}(1),1),...
                  size(local_answers{1}(2),2),...
                  sum(cellfun(@nnz,local_answers)));
        for n=1:numel(opspecs)
            A=A+local_answers{n};
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown operator type.');
        
end

end

% Input validation function
function grumble(spin_system,operators,spins)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~(ischar(operators)&&ischar(spins)))&&...
   (~(iscell(operators)&&iscell(spins)))&&...
   (~(ischar(operators)&&isnumeric(spins)))
    error('invalid operator specification.');
end
if iscell(operators)&&iscell(spins)&&(numel(operators)~=numel(spins))
    error('spins and operators cell arrays should have the same number of elements.');
end
if iscell(operators)&&any(~cellfun(@ischar,operators))
    error('all elements of the operators cell array should be strings.');
end
if iscell(spins)&&any(~cellfun(@isnumeric,spins))
    error('all elements of the spins cell array should be integer numbers.');
end
end

% It is nice to know that the computer understands the problem. But I would
% like to understand it too.
%
% Eugene Wigner

