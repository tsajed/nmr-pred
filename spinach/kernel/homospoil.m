% Emulates a strong homospoil pulse - only zero-frequency states
% survive the process. Syntax:
%
%                rho=homospoil(spin_system,rho,zqc_flag)
%
% Parameters:
%
%            rho - a state vector or a horizontal stack thereof
%
%       zqc_flag - a flag controlling the fate of zero-quantum 
%                  coherences. If set to 'keep', causes ZQCs to
%                  survive the process, approximating experimen-
%                  tal behaviour. If set to 'destroy', wipes the
%                  zero-quantum coherences - only the longitudi-
%                  nal states survive the process.
%
% Note: this function is only available for sphten-liouv formalism.
%
% i.kuprov@soton.ac.uk

function rho=homospoil(spin_system,rho,zqc_flag)

% Check consistency
grumble(spin_system,rho,zqc_flag)

% Store dimension statistics
spn_dim=size(spin_system.bas.basis,1);
spc_dim=numel(rho)/spn_dim;
problem_dims=size(rho);

% Fold indirect dimensions
rho=reshape(rho,[spn_dim spc_dim]);

% Pull the projection information from the basis
[~,M]=lin2lm(spin_system.bas.basis);

% Filter the state vector
switch zqc_flag
    
    case 'keep'
        
        % Find the states that have zero carrier frequency and kill everything else
        rho(abs(sum(repmat(spin_system.inter.basefrqs,size(spin_system.bas.basis,1),1).*M,2))>1e-6,:)=0;
    
    case 'destroy'
        
        % Find the longitudinal states and kill everything else
        rho(sum(abs(M),2)>0,:)=0;
    
    otherwise
        
        % Complain and bomb out
        error('unknown ZQC flag.');
        
end

% Unfold indirect dimensions
rho=reshape(rho,problem_dims);

% Report overly destructive calls
if norm(rho,1)<1e-10
    report(spin_system,'WARNING - all magnetization appears to have been destroyed by this call.');
end

end

% Consistency enforcement
function grumble(spin_system,rho,zqc_flag)
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available for sphten-liouv formalism.');
end
if ~isnumeric(rho)
    error('the state vector(s) must be numeric.');
end
if ~ischar(zqc_flag)
    error('zqc_flag parameter must be a character string.');
end
if ~ismember(zqc_flag,{'keep','destroy'})
    error('the available values for zqc_flag are ''keep'' and ''destroy''.');
end
end

% Rocket science has been mythologized all out of proportion to
% its true difficulty.
%
% John Carmack

