% Sets up interaction tensors under partial ordering in a liquid 
% crystal with the user-supplied order matrix. All adjustable pa-
% rameters are set during the call to create.m function. Syntax:
%
%               spin_system=residual(spin_system)
%
% Note: this function is only applicable to high-field NMR.
%
% Note: the function overwrites the interaction tensors supplied
%       by the user. Relaxation superoperator, if required, must
%       be computed before this function is called.
%
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk

function spin_system=residual(spin_system)

% Check consistency
grumble(spin_system);

% Process Zeeman interactions
for n=1:spin_system.comp.nspins
    if significant(spin_system.inter.zeeman.matrix{n},spin_system.tols.inter_cutoff)
        
        % Obtain isotropic part
        iso=eye(3)*trace(spin_system.inter.zeeman.matrix{n})/3;
        
        % Calculate residual order
        extra_zz=trace(spin_system.inter.order_matrix*(spin_system.inter.zeeman.matrix{n}-iso));
        
        % Update Zeeman tensor
        spin_system.inter.zeeman.matrix{n}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3]);
        
    end
end

% Process spin-spin couplings
for n=1:spin_system.comp.nspins
    for k=1:spin_system.comp.nspins
        if significant(spin_system.inter.coupling.matrix{n,k},spin_system.tols.inter_cutoff)
            
            % Obtain isotropic part
            iso=trace(spin_system.inter.coupling.matrix{n,k})*eye(3)/3;
            
            % Calculate residual order
            extra_zz=trace(spin_system.inter.order_matrix*(spin_system.inter.coupling.matrix{n,k}-iso));
            
            % Update coupling tensor
            spin_system.inter.coupling.matrix{n,k}=iso+diag([-extra_zz/3 -extra_zz/3 2*extra_zz/3]);
            
        end
    end
end

% Report back to the user
report(spin_system,'all interaction anisotropies have been replaced by their residuals.');

end

% Consistency enforcement
function grumble(spin_system)
if ~isfield(spin_system.inter,'order_matrix')
    error('order matrix infomation is missing from the spin_system structure.');
end
end

% If I'd observed all the rules, I'd never have got anywhere.
%
% Marilyn Monroe

