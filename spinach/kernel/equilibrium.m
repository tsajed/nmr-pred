% Returns the thermal equilibrium state at the current temperature. If the
% anisotropic part and the orientation parameters are not given, uses the
% isotropic Hamiltonian, otherwise uses the full Hamiltonian at the speci-
% fied orientation.
%
% If the temperature is set to zero during the call to create(), the high-
% temperature approximation to the thermal equilibrium state is returned. 
% If the temperature is specified, the accurate equilibrium state at that
% temperature is produced. Syntax:
%
%               rho=equilibrium(spin_system,H,Q,euler_angles)
%
% Arguments:
%
%    H            -  Isotropic part of the Hamiltonian left side pro-
%                    duct superoperator (in Lioville space) or Hamil-
%                    tonian (in Hilbert space).
%
%    Q            -  25 irreducible components of the anisotropic part
%                    of the Hamiltonian left side product superopera-
%                    tor (in Lioville space) or Hamiltonian (in Hil-
%                    bert space), as returned by hamiltonian.m
%
%    euler_angles -  a row vector of Euler angles (in radians) speci-
%                    fying the system orientation relative to the in-
%                    put orientation. If the angles are not supplied,
%                    only isotropic part of the Hamiltonian is used.
%
% WARNING: Liouville space calculations must supply left side product su-
%          peroperators, not commutation superoperators.
%
% WARNING: assumptions supplied to the hamiltonian.m call that generates
%          H and Q must be 'labframe'.
%
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk

function rho=equilibrium(spin_system,H,Q,euler_angles)

% Account for the orientation
if nargin==4
    
    % Check consistency
    grumble(H,Q,euler_angles);
    
    % Update Hamiltonian
    H=H+orientation(Q,euler_angles);
    
elseif nargin==2
    
    % Check consistency
    grumble(H);
    
else
    
    % Complain and bomb out
    error('incorrect number of input arguments.');
    
end

% Inform the user
report(spin_system,'computing the thermal equilibrium state...');

% Decide how to proceed
switch spin_system.bas.formalism
    
    case {'sphten-liouv','zeeman-liouv'}

        % Get the unit state
        unit=unit_state(spin_system);

        % Decide how to proceed
        if isempty(spin_system.rlx.temperature)
            
            % Inform the user
            report(spin_system,'WARNING: high-T approximation, density matrix trace ignored, normalization approximate.')
            
            % Use high-temperature approximation
            rho=-H*unit; rho=rho/norm(rho,1);
            
        else
            
            % Inform the user
            report(spin_system,['thermal equilibrium state at '...
                                num2str(spin_system.rlx.temperature) ' Kelvin.']);
            
            % Get temperature factor
            beta_factor=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature);
            
            % Get thermal equilibrium state
            rho=step(spin_system,H,unit,-1i*beta_factor); rho=rho/dot(unit,rho);
            
        end
        
    case {'zeeman-hilb'}
        
        % Decide how to proceed
        if isempty(spin_system.rlx.temperature)
            
            % Use high temperature approximation
            report(spin_system,'WARNING: high-T approximation, density matrix trace ignored, normalization approximate.');
            rho=-H; rho=rho/norm(rho,1);
            
        else
            
            % Inform the user
            report(spin_system,['thermal equilibrium state at '...
                                num2str(spin_system.rlx.temperature) ' Kelvin.']);
            
            % Get temperature factor
            beta_factor=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature);
            
            % Compute density matrix as a propagator in imaginary time
            rho=propagator(spin_system,H,-1i*beta_factor);
            
            % Match normalization to Liouville space
            rho=rho/trace(rho); rho=rho*sqrt(size(rho,1));
            
        end
        
end
    
end

% Consistency enforcement
function grumble(H,Q,euler_angles)
if ~isnumeric(H)
    error('Hamiltonian must be numeric.');
end
if nargin==3
    if (~iscell(Q))||any(size(Q)~=[5 5])||any(~cellfun(@isnumeric,Q(:)))
        error('Q must be a 5x5 cell array with numeric elements.');
    end
    if (~isnumeric(euler_angles))||(~isreal(euler_angles))||(numel(euler_angles)~=3)
        error('euler_angles must be a real three-element vector.');
    end
end
end

% Compassion is a wonderful thing. It's what one feels when one looks at
% a squashed caterpillar. An elevating experience. One can let oneself go
% and spread –- you know, like taking a girdle off. You don't have to hold
% your stomach, your heart or your spirit up –- when you feel compassion. 
% All you have to do is look down. It's much easier. When you look up, you
% get a pain in the neck. [...] Oh, it has an antithesis – but such a hard,
% demanding one... admiration, Mrs. Jones, admiration. But that takes more
% than a girdle. So I say that anyone for whom we can't feel sorry is a
% vicious person.
%
% Ayn Rand, "The Fountainhead"

