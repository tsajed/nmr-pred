% Basis set control. This is the second mandatory function (after create.m)
% that must be executed in every calculation to get the Spinach kernel go-
% ing. See the Basis Selection section of the Spinach manual for the detail-
% ed description of the various options. 
%
% Note: it is very important that you understand the factors that influence
% basis set selection in spin dynamics simulations -- see the following pa-
% per for further details on the subject:
%
%             http://link.aip.org/link/doi/10.1063/1.3624564
%
% i.kuprov@soton.ac.uk

function spin_system=basis(spin_system,bas)

% Show the banner
banner(spin_system,'basis_banner');

% Check the input
grumble(spin_system,bas);

% Report the basis to the user
spin_system.bas.formalism=bas.formalism;
switch spin_system.bas.formalism
    case 'zeeman-hilb'
        report(spin_system,'Zeeman basis set using Hilbert space matrix formalism.');
    case 'zeeman-liouv'
        report(spin_system,'Zeeman basis set using Liouville space matrix formalism.');
    case 'sphten-liouv'
        report(spin_system,'spherical tensor basis set using Liouville space matrix formalism.');
    otherwise
        error('unrecognized formalism - see the basis preparation section of the manual.');
end

% Save bas.approximation field to the spin_system
spin_system.bas.approximation=bas.approximation;

% Process spherical tensor basis sets
if strcmp(spin_system.bas.formalism,'sphten-liouv')

    % Report approximation level to the user
    switch spin_system.bas.approximation
        case 'IK-0'
            report(spin_system,['IK-0 approximation - all correlations of all spins up to order ' num2str(bas.level)]);
        case 'IK-1'
            report(spin_system,['IK-1 approximation - spin correlations up to order ' num2str(bas.level) ' between directly coupled spins.']);
            report(spin_system,['IK-1 approximation - spin correlations up to order ' num2str(bas.space_level) ' between all spins within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
        case 'IK-2'
            report(spin_system, 'IK-2 approximation - spin correlations involving all nearest neighbours of each spin on the coupling graph.');
            report(spin_system,['IK-2 approximation - spin correlations up to order ' num2str(bas.space_level) ' between all spins within ' num2str(spin_system.tols.prox_cutoff) ' Angstrom of each other.']);
        case 'ESR-1'
            report(spin_system, 'ESR-1 approximation - complete basis on all electrons.');
            report(spin_system, 'ESR-1 approximation - zero-quantum basis on all nuclei.');
        case 'ESR-2'
            report(spin_system, 'ESR-2 approximation - complete basis on all electrons.');
            report(spin_system, 'ESR-2 approximation - complete basis on all nuclei with anisotropic HFC tensors.');
            report(spin_system, 'ESR-2 approximation - zero-quantum basis on all nuclei with isotropic HFC tensors.');
        case 'none'
            report(spin_system, 'complete basis set on all spins.');
        otherwise
            error('unrecognized approximation level - see the basis set preparation manual.');
    end

    % Run connectivity analysis
    if ismember(spin_system.bas.approximation,{'IK-1','IK-2'})
        
        % Build connectivity and proximity matrices
        switch bas.connectivity
            case 'scalar_couplings'
                report(spin_system,'scalar couplings will be used to build the coupling graph.');
                spin_system.inter.conmatrix=sparse(abs(cellfun(@trace,spin_system.inter.coupling.matrix)/3)>spin_system.tols.inter_cutoff);
            case 'full_tensors'
                report(spin_system,'full coupling tensors will be used to build the coupling graph.');
                spin_system.inter.conmatrix=sparse(cellfun(@norm,spin_system.inter.coupling.matrix)>spin_system.tols.inter_cutoff);
        end
        report(spin_system,['coupling tensors with norm below ' num2str(spin_system.tols.inter_cutoff) ' rad/s will be ignored.']);
        
        % Make sure each spin is connected and proximate to itself
        spin_system.inter.conmatrix=spin_system.inter.conmatrix|speye(size(spin_system.inter.conmatrix));
        spin_system.inter.proxmatrix=spin_system.inter.proxmatrix|speye(size(spin_system.inter.proxmatrix));
        
        % Make sure connectivity and proximity are reciprocal
        spin_system.inter.conmatrix=spin_system.inter.conmatrix|transpose(spin_system.inter.conmatrix);
        spin_system.inter.proxmatrix=spin_system.inter.proxmatrix|transpose(spin_system.inter.proxmatrix);
        
        % Issue a report to the user
        report(spin_system,['connectivity matrix density ' num2str(100*nnz(spin_system.inter.conmatrix)/numel(spin_system.inter.conmatrix)) '%']);
        report(spin_system,['proximity matrix density ' num2str(100*nnz(spin_system.inter.proxmatrix)/numel(spin_system.inter.proxmatrix)) '%']);
        
        % Determine the number of independent subsystems
        n_subsystems=max(scomponents(spin_system.inter.conmatrix|spin_system.inter.proxmatrix));
        
        % Print a notice to the user
        if n_subsystems>1
            report(spin_system,['WARNING - the system contains ' num2str(n_subsystems) ' non-interacting subsystems.']);
        end
        
    end

    % Generate subgraphs for IK-0 basis
    if strcmp(spin_system.bas.approximation,'IK-0')
        
        % Find all possible groups of bas.level spins
        col_index=combnk(1:spin_system.comp.nspins,bas.level);
        
        % Get the number of groups
        ngroups=size(col_index,1);
        
        % Assign numbers to groups
        row_index=kron((1:ngroups)',ones(1,bas.level));
        
        % Generate the subgraph list
        coupling_subgraphs=logical(sparse(row_index,col_index,ones(size(col_index)),ngroups,spin_system.comp.nspins));
        report(spin_system,['' num2str(size(coupling_subgraphs,1)) ' subgraphs generated by combinatorial analysis.']);
        
        % Do not run proximity analysis
        proximity_subgraphs=[];
        
    end
    
    % Generate subgraphs for IK-1 basis
    if strcmp(spin_system.bas.approximation,'IK-1')
        
        % Run connectivity analysis
        coupling_subgraphs=dfpt(spin_system.inter.conmatrix,bas.level);
        report(spin_system,['' num2str(size(coupling_subgraphs,1)) ' subgraphs generated from coupling data.']);
        
        % Run proximity analysis
        proximity_subgraphs=dfpt(spin_system.inter.proxmatrix,bas.space_level);
        report(spin_system,['' num2str(size(proximity_subgraphs,1)) ' subgraphs generated from proximity data.']);
        
    end
    
    % Generate subgraphs for IK-2 basis
    if strcmp(spin_system.bas.approximation,'IK-2')
        
        % Run connectivity analysis
        coupling_subgraphs=spin_system.inter.conmatrix;
        report(spin_system,['' num2str(size(coupling_subgraphs,1)) ' subgraphs generated from coupling data.']);
        
        % Run proximity analysis
        proximity_subgraphs=dfpt(spin_system.inter.proxmatrix,bas.space_level);
        report(spin_system,['' num2str(size(proximity_subgraphs,1)) ' subgraphs generated from proximity data.']);
        
    end

    % Build IK-0,1,2 basis sets
    if ismember(spin_system.bas.approximation,{'IK-0','IK-1','IK-2'})
        
        % Include user-specified subgraphs
        if (~isfield(bas,'manual'))||isempty(bas.manual)
            manual_subgraphs=[];
        else
            manual_subgraphs=bas.manual;
            report(spin_system,['added ' num2str(size(manual_subgraphs,1)) ' subgraphs specified by the user.']);
        end
        
        % Assemble the subgraph list
        subgraphs=[coupling_subgraphs; proximity_subgraphs; manual_subgraphs];
        clear('coupling_subgraphs','proximity_subgraphs','manual_subgraphs');

        % Prune subgraphs involving spin zero particles
        report(spin_system,'pruning subgraphs involving zero spin particles...');
        subgraphs(:,spin_system.comp.mults==1)=0;
        
        % Remove identical subgraphs
        report(spin_system,'removing identical subgraphs...');
        subgraphs=unique(subgraphs,'rows');

        % Store subgraphs for future use
        spin_system.bas.subgraphs=logical(subgraphs);
        
        % Report back to the user
        subgraph_sizes=sum(subgraphs,2);
        for n=min(subgraph_sizes):max(subgraph_sizes)
            if nnz(subgraph_sizes==n)>0
                report(spin_system,['generated ' num2str(nnz(subgraph_sizes==n)) ' subgraphs of size ' num2str(n)]);
            end
        end
        
        % Balance the subgraph list
        subgraphs=subgraphs(randperm(size(subgraphs,1)),:);
        
        % Populate the basis descriptor array
        report(spin_system,'building basis descriptor...'); 
        spmd
            
            % Localize the problem at the nodes
            subgraphs=codistributed(subgraphs,codistributor('1d',1));
            subgraphs_loc=getLocalPart(subgraphs);
        
            % Preallocate local basis descriptor array
            nstates=0; 
            for n=1:size(subgraphs_loc,1)
                nstates=nstates+prod(spin_system.comp.mults(logical(subgraphs_loc(n,:))))^2;
            end
            basis_spec=spalloc(nstates,spin_system.comp.nspins,nstates*max(subgraph_sizes));
        
            % Populate local basis descriptor array
            list_position=1;
            for n=1:size(subgraphs_loc,1)
                
                % Determine the total number of states in the current subgraph
                nstates=prod(spin_system.comp.mults(logical(subgraphs_loc(n,:))))^2;
                
                % Determine which spins belong to the current subgraph
                spins_involved=find(subgraphs_loc(n,:));
                
                % Build the basis table
                for k=1:numel(spins_involved)
                    basis_spec(list_position:(list_position+nstates-1),spins_involved(k))=...
                    kron(kron(ones(prod(spin_system.comp.mults(spins_involved(1:(k-1))))^2,1),...
                    (0:(spin_system.comp.mults(spins_involved(k))^2-1))'),...
                    ones(prod(spin_system.comp.mults(spins_involved((k+1):end)))^2,1));
                end
                list_position=list_position+nstates;
                
            end
            
            % Prune local basis descriptor array
            basis_spec=unique(basis_spec,'rows');
            
        end
        
        % Pull the basis descriptor from the nodes
        try
            basis_spec=vertcat(basis_spec{:});
        catch
            error('parpool appears to be stuck, check your parallel configuration.');
        end
        
    end

    % Build ESR-1 basis set
    if strcmp(spin_system.bas.approximation,'ESR-1')
        
        % Prepare direct product components
        components=cell(1,spin_system.comp.nspins);
        for n=1:spin_system.comp.nspins
            if strcmp(spin_system.comp.isotopes{n}(1),'E')
                components{n}=(0:(spin_system.comp.mults(n)^2-1))';
                report(spin_system,['ESR-1 approximation, complete basis set on spin ' num2str(n)]);
            else
                components{n}=((1:spin_system.comp.mults(n)).*((1:spin_system.comp.mults(n))-1))';
                report(spin_system,['ESR-1 approximation, zero-quantum basis set on spin ' num2str(n)]);
            end
        end
        
        % Preallocate basis table
        basis_spec=zeros(prod(cellfun(@numel,components)),spin_system.comp.nspins);
        
        % Generate basis specification
        report(spin_system,'building basis descriptor...');
        parfor n=1:spin_system.comp.nspins
            dimension_before=prod(cellfun(@numel,components(1:(n-1)))); %#ok<PFBNS>
            dimension_after=prod(cellfun(@numel,components((n+1):end)));
            basis_spec(:,n)=kron(kron(ones(dimension_before,1),...
                            components{n}),ones(dimension_after,1));
        end
        
    end

    % Process ESR-2 basis set
    if strcmp(spin_system.bas.approximation,'ESR-2')
        
        % Analyze coupling structure
        spins_with_full_basis=false(1,spin_system.comp.nspins);
        for n=1:spin_system.comp.nspins
            
            % Use complete basis set on anything that has an anisotropic coupling
            for k=1:spin_system.comp.nspins
                if (norm(spin_system.inter.coupling.matrix{n,k})>spin_system.tols.inter_cutoff)&&...
                    norm(spin_system.inter.coupling.matrix{n,k}-eye(3)*trace(spin_system.inter.coupling.matrix{n,k})/3)>spin_system.tols.inter_cutoff
                    spins_with_full_basis(n)=true(); spins_with_full_basis(k)=true();
                end
            end
            
            % Use complete basis set on all electrons
            if strcmp(spin_system.comp.isotopes{n}(1),'E')
                spins_with_full_basis(n)=true();
            end
            
        end
        
        % Prepare direct product components
        components=cell(1,spin_system.comp.nspins);
        for n=1:spin_system.comp.nspins
            if spins_with_full_basis(n)
                components{n}=(0:(spin_system.comp.mults(n)^2-1))';
                report(spin_system,['ESR-2 approximation - complete basis set on spin ' num2str(n)]);
            else
                components{n}=((1:spin_system.comp.mults(n)).*((1:spin_system.comp.mults(n))-1))';
                report(spin_system,['ESR-2 approximation - zero-quantum basis set on spin ' num2str(n)]);
            end
        end
        
        % Preallocate the basis table
        basis_spec=zeros(prod(cellfun(@numel,components)),spin_system.comp.nspins);
        
        % Generate the basis specification
        report(spin_system,'building basis descriptor...');
        parfor n=1:spin_system.comp.nspins
            dimension_before=prod(cellfun(@numel,components(1:(n-1)))); %#ok<PFBNS>
            dimension_after=prod(cellfun(@numel,components((n+1):end)));
            basis_spec(:,n)=kron(kron(ones(dimension_before,1),...
                            components{n}),ones(dimension_after,1));
        end
        
    end

    % Process complete basis set
    if strcmp(spin_system.bas.approximation,'none')
        
        % Preallocate basis table
        basis_spec=zeros(prod(spin_system.comp.mults)^2,spin_system.comp.nspins);
        
        % Build basis table
        report(spin_system,'building basis descriptor...');
        parfor k=1:spin_system.comp.nspins
            basis_spec(:,k)=kron(kron(ones(prod(spin_system.comp.mults(1:(k-1)))^2,1),...
                            (0:(spin_system.comp.mults(k)^2-1))'),...
                            ones(prod(spin_system.comp.mults((k+1):end))^2,1)); %#ok<PFBNS>
        end
        
    end

    % Apply coherence order filter
    if isfield(bas,'projections')
        
        % Inform the user
        report(spin_system,['applying coherence order filter with M=[' num2str(bas.projections) ']...']);
        
        % Compute coherence order for each basis element
        [~,M]=lin2lm(basis_spec);
        
        % Keep specified coherence orders and the unit state
        state_mask=false(size(basis_spec,1),1); 
        state_mask(1)=true(); projection_numbers=sum(M,2);
        for n=bas.projections
            state_mask=state_mask|(projection_numbers==n);
        end
        basis_spec(~state_mask,:)=[];
        
    end
    
    % Apply zero-quantum filter
    if isfield(bas,'longitudinals')
        
        % Start with empty mask
        state_mask=false(size(basis_spec,1),1);
        
        % Loop over the specified spins
        for n=1:numel(bas.longitudinals)
            
            % Inform the user
            report(spin_system,['applying longitudinal spin order filter on ' bas.longitudinals{n}]);
            
            % Kill the specified transverse states
            for k=find(strcmp(bas.longitudinals{n},spin_system.comp.isotopes))
                
                % Analyze the basis
                [~,M]=lin2lm(basis_spec(:,k));
                
                % Update the mask
                state_mask=or(state_mask,M~=0);
                
            end
            
        end
        
        % Kill the states
        basis_spec(state_mask,:)=[];
        
    end
    
    % Remove spin correlations between chemical subsystems
    report(spin_system,'removing correlations between chemical reaction endpoints...');
    for n=1:numel(spin_system.chem.parts)
        for k=1:numel(spin_system.chem.parts)
            if n~=k
                
                % Find states involving the first subsystem
                spins_a=spin_system.chem.parts{n};
                states_a=logical(sum(basis_spec(:,spins_a),2));
                
                % Find states involving the second subsystem
                spins_b=spin_system.chem.parts{k};
                states_b=logical(sum(basis_spec(:,spins_b),2));
                
                % Remove the intersection
                basis_spec(states_a&states_b,:)=[];
                
            end
        end
    end

    % Write the final basis specification
    report(spin_system,'sorting the basis...');
    spin_system.bas.basis=unique(basis_spec,'rows');
    
    % Build state-cluster cross-membership list if needed
    if ismember('xmemlist',spin_system.sys.enable)
        
        % Report to the user
        report(spin_system,'building state-subgraph cross-membership list... ');
        
        % Localise variables
        basis_loc=spin_system.bas.basis;
        subgraphs_loc=spin_system.bas.subgraphs;
        
        % Preallocate the result
        nstates=size(basis_loc,1); nclusters=size(subgraphs_loc,1);
        xmemlist=spalloc(nstates,nclusters,ceil(nstates*nclusters/spin_system.comp.nspins));
        
        % Run the matching
        parfor n=1:size(subgraphs_loc,1)
            xmemlist(:,n)=~any(basis_loc(:,~subgraphs_loc(n,:)),2); %#ok<SPRIX,PFBNS>
        end
        
        % Store the result
        spin_system.bas.xmemlist=xmemlist;

    end 
    
    % Print the summary
    summary(spin_system,'basis');

    % Run the symmetry treatment
    spin_system=symmetry(spin_system,bas);

end

% Process Hilbert space Zeeman basis
if strcmp(spin_system.bas.formalism,'zeeman-hilb')
   
    % Preallocate basis set array
    spin_system.bas.basis=zeros(prod(spin_system.comp.mults),spin_system.comp.nspins);
    
    % Fill basis set array
    for n=1:spin_system.comp.nspins
        current_column=1;
        for k=1:spin_system.comp.nspins
            if n==k
                current_column=kron(current_column,(1:spin_system.comp.mults(k))');
            else
                current_column=kron(current_column,ones(spin_system.comp.mults(k),1));
            end
        end
        spin_system.bas.basis(:,n)=current_column;
    end 
    
    % Run the symmetry treatment
    spin_system=symmetry(spin_system,bas);
    
end

end

% Grumble function
function grumble(spin_system,bas)

% Check bas.formalism
if ~isfield(bas,'formalism')
    error('basis specification in bas.formalism is required.');
elseif ~ischar(bas.formalism)
    error('bas.formalism must be a string.');
elseif ~ismember(bas.formalism,{'zeeman-hilb','zeeman-liouv','sphten-liouv'})
    error('unrecognized formalism - see the basis preparation section of the manual.');
end

% Check zeeman-hilb formalism options
if strcmp(bas.formalism,'zeeman-hilb')
    
    % Check bas.approximation
    if ~isfield(bas,'approximation')
        error('approximation level must be specified in bas.approximation for zeeman-hilb formalism.');
    elseif ~ischar(bas.approximation)
        error('bas.approximation must be a string.');
    elseif ~ismember(bas.approximation,{'none'})
        error('bas.approximation should be set to ''none'' in zeeman-hilb formalism.');
    end

end

% Check zeeman-liouv formalism options
if strcmp(bas.formalism,'zeeman-liouv')
    
    % Check bas.approximation
    if ~isfield(bas,'approximation')
        error('approximation level must be specified in bas.approximation for zeeman-liouv formalism.');
    elseif ~ischar(bas.approximation)
        error('bas.approximation must be a string.');
    elseif ~ismember(bas.approximation,{'none'})
        error('bas.approximation should be set to ''none'' in zeeman-liouv formalism.');
    end

end

% Check sphten-liouv formalism options
if strcmp(bas.formalism,'sphten-liouv')
    
    % Check bas.approximation
    if ~isfield(bas,'approximation')
        error('approximation level must be specified in bas.approximation for sphten-liouv formalism.');
    elseif ~ischar(bas.approximation)
        error('bas.approximation must be a string.');
    elseif ~ismember(bas.approximation,{'IK-0','IK-1','IK-2','ESR-1','ESR-2','none'})
        error('unrecognized approximation - see the basis preparation section of the manual.');
    end
    
    % Check bas.connectivity
    if ismember(bas.approximation,{'IK-1','IK-2'})
        if ~isfield(bas,'connectivity')
            error('connectivity type must be specified in bas.connectivity variable.');
        elseif ~ischar(bas.connectivity)
            error('bas.connectivity must be a string.');
        elseif ~ismember(bas.connectivity,{'scalar_couplings','full_tensors'})
            error('unknown connectivity type - see the basis preparation section of the manual.');
        end
    end
    
    % Check bas.level
    if ismember(bas.approximation,{'IK-0','IK-1'})&&(~isfield(bas,'level'))
        error('connectivity tracing depth must be specified in bas.level variable.');
    end
    if isfield(bas,'level')
        if (~isnumeric(bas.level))||(~isscalar(bas.level))||(mod(bas.level,1)~=0)||(bas.level<1)
            error('bas.level must be a positive integer.');
        end
    end
    
    % Check bas.space_level
    if ismember(bas.approximation,{'IK-1','IK-2'})&&(~isfield(bas,'space_level'))
        error('proximity tracing depth must be specified in bas.space_level variable.');
    end
    if isfield(bas,'space_level')
        if  (~isnumeric(bas.space_level))||(~isscalar(bas.space_level))||(mod(bas.space_level,1)~=0)||(bas.space_level<1)
            error('bas.space_level must be a positive integer.');
        end
    end
    
    % Check bas.manual
    if isfield(bas,'manual')
        if (~islogical(bas.manual))&&(~isnumeric(bas.manual))
            error('bas.manual must be a logical matrix.');
        elseif size(bas.manual,2)~=spin_system.comp.nspins
            error('the number of rows in bas.manual must be equal to the number of spins in the system.')
        end
    end
    
    % Check bas.projections
    if isfield(bas,'projections')
        if (~isnumeric(bas.projections))||(~isrow(bas.projections))||any(mod(bas.projections,1)~=0)
            error('bas.projections must be a row vector of integers.');
        end
    end
    
    % Check bas.longitudinals
    if isfield(bas,'longitudinals')
        if (~iscell(bas.longitudinals))||any(~cellfun(@ischar,bas.longitudinals))
            error('bas.longitudinals must be a cell array of strings.');
        end
        if any(~ismember(bas.longitudinals,spin_system.comp.isotopes))
            error('bas.longitudinals refers to spins that are not present in the system.');
        end
    end
    
end
    
% Disallow approximation options with exact formalism
if ~ismember(bas.approximation,{'IK-0','IK-1','IK-2'})
    if isfield(bas,'level')||isfield(bas,'space_level')
        error('bas.level and bas.space_level are only applicable to IK-0,1,2 basis sets.');
    end
end

% Enforce sphten-liouv with projection selection
if isfield(bas,'projections')&&(~strcmp(bas.formalism,'sphten-liouv'))
    error('bas.projections option is only available for sphten-liouv formalism.');
end

% Enforce sphten-liouv with chemical exchange
if (~isscalar(spin_system.chem.rates))&&(~strcmp(bas.formalism,'sphten-liouv'))
    error('chemical reaction modelling is only available for sphten-liouv formalism.');
end

end

% In 1969, Robert Rathbun Wilson, the US physicist who headed Fermilab, the world's
% highest-energy particle accelerator laboratory, addressed the Congressional Joint
% Committee on Atomic Energy. Rhode Island Senator John Pastore asked Wilson to spell
% out what research into high-energy particle physics would do to improve the defence
% of the United States. Wilson gave a reply that went down in scientific history. Fer-
% milab, he said, had "nothing to do directly with defending our country, except to
% make it worth defending".
%
% http://www.theregister.co.uk/2009/02/09/woudhuysen_energise_1/

