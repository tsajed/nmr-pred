% Modifies a symmetric relaxation superoperator to relax the system
% towards a user-specified state. Liouville space spherical tensor
% formalism only. Syntax:
%
%                  R=thermalize(spin_system,R,rho)
%
% i.kuprov@soton.ac.uk

function R=thermalize(spin_system,R,rho)

if strcmp(spin_system.bas.formalism,'sphten-liouv')
    
    % Modify the relaxation superoperator
    R(1,1)=1; R(:,1)=-R*rho;
    
else
    
    % Complain and bomb out
    error('this function is only available in sphten-liov formalism.');
    
end

end

% As for the review experience, remember that story about 
% Pavlov's dogs? Conditional reflexes and things. There's
% a less well known experiment when dogs are punished and
% rewarded in a way that's uncorrelated with what they do.
% The dogs eventually develop schizophrenia... that's how
% that review made me feel.
%
% (from IK's email to the project team,
%  after the final review meeting on an
%  EU grant, in which the outcomes were
%  praised by Brussels to high heaven)

