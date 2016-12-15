% Liquid-phase interface to pulse sequences. Generates a Liouvillian
% superoperator and passes it on to the pulse sequence function, which
% should be supplied as a handle.
%
% This interface handles RDC mode -- if parameters.rdc switch is set,
% it would use the order matrix supplied by the user to compute the
% residual anisotropies of all interactions. Syntax:
%
%    fid=liquid(spin_system,pulse_sequence,parameters,assumptions)
%
% Arguments:
%
%    pulse_sequence  -  pulse sequence function handle. See the
%                       experiments directory for the list of
%                       pulse sequences that ship with Spinach.
% 
%    parameters      -  a structure containing sequence-specific
%                       parameters. See the pulse sequence header
%                       for the list of parameters each sequence
%                       requires. The parameters that this inter-
%                       face itself requires are:
%                      
%                       parameters.spins - a cell array giving the
%                       spins that the pulse sequence works on, in
%                       the order of channels, e.g. {'1H','13C'}
%
%                       parameters.offset - a cell array giving
%                       transmitter offsets on each of the spins
%                       listed in parameters.spins array.
%
%                       parameters.rdc - a logical switch control-
%                       ling the inclusion of residual dipolar co-
%                       uplings into the calculation.
%
%    assumptions     -  context-specific assumptions ('nmr', 'epr',
%                       'labframe', etc.) - see the pulse sequence
%                       header for information on this setting.
%
% This function returns whatever it is that the pulse sequence returns.
%
% i.kuprov@soton.ac.uk

function answer=liquid(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Set assumptions
spin_system=assume(spin_system,assumptions);

% Get the Liouvillian
if isfield(parameters,'rdc')&&parameters.rdc
    
    % With RDCs, first get the relaxation and kinetics
    R=relaxation(spin_system); K=kinetics(spin_system);
    
    % Then do the liquid crystal averaging
    spin_system=residual(spin_system);
    
    % Then get the coherent parts of the Liouvillian
    [I,Q]=hamiltonian(spin_system); H=I+orientation(Q,[0 0 0]);

else
    
    % Get the isotropic Liouvillian
    H=hamiltonian(spin_system);
    R=relaxation(spin_system);
    K=kinetics(spin_system);

end

% Process channel offsets
H=frqoffset(spin_system,H,parameters);

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Apply rotating frames
for k=1:numel(parameters.rframes)
    H=rotframe(spin_system,C{k},H,parameters.rframes{k}{1},parameters.rframes{k}{2});
end

% Get problem dimensions
parameters.spc_dim=1; parameters.spn_dim=size(H,1);

% Report to the user
report(spin_system,'running the pulse sequence...');

% Call the pulse sequence
answer=pulse_sequence(spin_system,parameters,H,R,K);

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'decouple')
    report(spin_system,'parameters.decouple field not set, assuming no decoupling.');
    parameters.decouple={};
end
if ~isfield(parameters,'rframes')
    report(spin_system,'parameters.rframes field not set, assuming no additional rotating frame transformations.');
    parameters.rframes={};
end
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'verbose')
    report(spin_system,'parameters.verbose field not set, silencing array operations.');
    parameters.verbose=0;
end
end

% Consistency enforcement
function grumble(spin_system,pulse_sequence,parameters,assumptions)

% Pulse sequence
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end

% Assumptions
if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
end

% Active spins
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end

% Offsets
if isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins variable must be specified alongside parameters.offset variable.');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset variable must have the same number of elements as parameters.spins.');
end

% Rotating frame transformations
if ~isfield(parameters,'rframes')
    error('parameters.rframes variable must be specified.');
elseif ~iscell(parameters.rframes)
    error('parameters.rframes must be a cell array.');
end
for n=1:numel(parameters.rframes)
    if ~iscell(parameters.rframes{n})
        error('elements of parameters.rframes must be cell arrays.');
    end
    if numel(parameters.rframes{n})~=2
        error('elements of parameters.rframes must have exactly two sub-elements each.');
    end
    if ~ischar(parameters.rframes{n}{1})
        error('the first part of each element of parameters.rframes must be a character string.');
    end
    if ~ismember(parameters.rframes{n}{1},spin_system.comp.isotopes)
        error('parameters.rframes refers to a spin that is not present in the system.');
    end
    if ~isnumeric(parameters.rframes{n}{2})
        error('the second part of each element of parameters.rframes must be a number.');
    end
end

end

% A man who can conceive a thing as beautiful as this should never 
% have allowed it to be erected. [...] But he will let it be built,
% so that women will hang out diapers on his terraces, so that men
% will spit on his stairways and draw dirty pictures on his walls.
% [...] They will be committing only a mean little indecency, but 
% he has committed a sacrilege.
%
% Ayn Rand, "The Fountainhead"

