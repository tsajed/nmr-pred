% Powder interface to pulse sequences. Generates a Liouvillian super-
% operator, the initial state and the coil state, then passes them on
% to the pulse sequence function.
%
% This function supports parallel processing via Matlab's Distributed
% Computing Toolbox - different crystal orientations are evaluated on
% different labs. Arguments:
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
%                       listed in parameters.spins.
%
%                       parameters.grid - name of the spherical 
%                       averaging grid file (see the grids direc-
%                       tory in the kernel).
%
%    assumptions     -  context-specific assumptions ('nmr', 'epr',
%                       'labframe', etc.) - see the pulse sequence
%                       header for information on this setting.
%
% This function returns a powder average of whatever it is that the pul-
% se sequence returns. If a structure is returned by the pulse sequence,
% the structures are powder averaged field-by-field.
%
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk

function answer=powder(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Set the assumptions
spin_system=assume(spin_system,assumptions);

% Get the rotational expansion of the Hamiltonian
[I,Q]=hamiltonian(spin_system);

% Get kinetics superoperator
K=kinetics(spin_system);

% Add offsets to the isotropic part
I=frqoffset(spin_system,I,parameters);

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Get problem dimensions
parameters.spc_dim=1; parameters.spn_dim=size(I,1);

% Get the averaging grid as a struct
spherical_grid=load([spin_system.sys.root_dir '/kernel/grids/' parameters.grid],...
                                           'alphas','betas','gammas','weights');
% Assign local variables
alphas=spherical_grid.alphas; betas=spherical_grid.betas; 
gammas=spherical_grid.gammas; weights=spherical_grid.weights;

% Preallocate answer array
ans_array=cell(numel(weights),1);

% Shut up and inform the user
report(spin_system,['powder average being computed over ' num2str(numel(weights)) ' orientations.']);
if ~isfield(parameters,'verbose')||(parameters.verbose==0)
    report(spin_system,'pulse sequence has been silenced to avoid excessive output.')
    spin_system.sys.output='hush';
end

% Crunch orientations in parallel
parfor n=1:numel(weights)
    
    % Get the full Hamiltonian at the current orientation
    H=I+orientation(Q,[gammas(n) betas(n) alphas(n)]); H=(H+H')/2;
    
    % Apply rotating frames
    for k=1:numel(parameters.rframes) %#ok<PFBNS>
        H=rotframe(spin_system,C{k},H,parameters.rframes{k}{1},parameters.rframes{k}{2}); %#ok<PFBNS>
    end
    
    % Get the relaxation superoperator at the current orientation
    R=relaxation(spin_system,[gammas(n) betas(n) alphas(n)]);
    
    % Report to the user
    report(spin_system,'running the pulse sequence...');
    
    % Run the simulation (it might return a structure)
    ans_array{n}=pulse_sequence(spin_system,parameters,H,R,K); %#ok<PFBNS>
    
end

% Add up structures
answer=weights(1)*ans_array{1};
for n=2:numel(ans_array)
    answer=answer+weights(n)*ans_array{n};
end

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

% Consistency checking
function grumble(spin_system,pulse_sequence,parameters,assumptions)

% Spherical grid
if ~isfield(parameters,'grid')
    error('spherical averaging grid must be specified in parameters.grid variable.');
elseif isempty(parameters.grid)
    error('parameters.grid variable cannot be empty.');
elseif ~ischar(parameters.grid)
    error('parameters.grid variable must be a character string.');
end

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

% Do not confuse altruism with kindness. The irreducible primary of
% altruism is self-sacrifice - which means self-immolation, self-
% abnegation, self-denial. Do not hide behind such superficialities
% as whether you should or should not give a dime to a beggar. This
% is not the issue. The issue is whether you do or do not have the
% right to exist without giving him that dime. The issue is whether
% the need of others is the first mortgage on your life and the mo-
% ral purpose of your existence. Any man of self-esteem will answer:
% No. Altruism says: Yes.
%
% Ayn Rand

