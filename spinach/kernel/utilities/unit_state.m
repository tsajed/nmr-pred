% Returns a unit state in the current formalism and basis. Syntax:
%
%                    rho=unit_state(spin_system)
%
% There are no adjustable parameters.
%
% i.kuprov@soton.ac.uk
% d.savostyanov@soton.ac.uk

function rho=unit_state(spin_system)

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Normalized T(0,0) state
        rho=sparse(1,1,1,size(spin_system.bas.basis,1),1);
        
    case 'zeeman-liouv'
        
        % Normalized stretched unit matrix
        rho=speye(prod(spin_system.comp.mults));
        rho=rho(:); rho=rho/norm(rho,2);
        
    case 'zeeman-hilb'
        
        % Unit matrix
        rho=speye(prod(spin_system.comp.mults));
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% There used to be a simple story about Russian literature, that we
% thought the good writers were the ones who opposed the regime. Once
% we don't have that story about Russia as a competitor, or an enemy,
% it was much less clear to us what we should be interested in.
%
% Edwin Frank, the editor of NYRB Classics

