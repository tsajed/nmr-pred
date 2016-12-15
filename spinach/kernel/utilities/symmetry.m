% Symmetry treatment. This is a service function of the Spinach
% kernel that should not be called directly.
%
% i.kuprov@soton.ac.uk
% hannah.hogben@chem.ox.ac.uk

function spin_system=symmetry(spin_system,bas)

% Check consistency
grumble(spin_system,bas);

% Check the disable switch
if ismember('symmetry',spin_system.sys.disable)
    
    % Issue a reminder to the user
    report(spin_system,'WARNING - symmetry factorization disabled by the user.');
    
    % Write empty cells
    spin_system.comp.sym_group={};
    spin_system.comp.sym_spins={};
    spin_system.comp.sym_a1g_only=true();
    
else

    % Symmetry group
    if isfield(bas,'sym_group')
        spin_system.comp.sym_group=bas.sym_group;
    else
        spin_system.comp.sym_group={};
    end
    
    % Symmetry-related spins
    if isfield(bas,'sym_spins')
        spin_system.comp.sym_spins=bas.sym_spins;
    else
        spin_system.comp.sym_spins={};
    end
    
    % Irreducible representation composition
    if isfield(bas,'sym_a1g_only')
        spin_system.comp.sym_a1g_only=bas.sym_a1g_only;
    else
        spin_system.comp.sym_a1g_only=true();
    end
    
    % Report back to the user
    if ~isempty(spin_system.comp.sym_group)
        summary(spin_system,'symmetry','permutation symmetry summary');
    end
    
    % Compute group direct product if necessary
    if numel(spin_system.comp.sym_group)>1
        
        % Lift constituent groups from the database
        ngroups=numel(spin_system.comp.sym_group); groups=cell(1,ngroups);
        for n=1:ngroups
            groups{n}=perm_group(spin_system.comp.sym_group{n});
        end
        
        % Compute direct product character table
        group.characters=1;
        for n=1:ngroups
            group.characters=kron(group.characters,groups{n}.characters);
        end
        group.irrep_dims=group.characters(:,1)';
        group.n_irreps=size(group.characters,1);
        report(spin_system,['' num2str(group.n_irreps) ' irreps in the group direct product.']);
        report(spin_system,['dimensions of the irreps ' num2str(group.irrep_dims)]);
        
        % Compute direct product element list
        group.elements=groups{1}.elements; group.order=groups{1}.order;
        for n=2:ngroups
            group.elements=[kron(group.elements,ones(groups{n}.order,1))...
                kron(ones(group.order,1),groups{n}.elements+size(group.elements,2))];
            group.order=group.order*groups{n}.order;
        end
        report(spin_system,['' num2str(size(group.elements,1)) ' symmetry operations in the group direct product.']);
        
        % Concatenate spin lists
        spins=horzcat(spin_system.comp.sym_spins{:});
        
    elseif numel(spin_system.comp.sym_group)==1
        
        % Lift the group from the database
        spins=spin_system.comp.sym_spins{1};
        group=perm_group(spin_system.comp.sym_group{1});
        report(spin_system,['' num2str(group.n_irreps) ' irreps in the symmetry group.']);
        report(spin_system,['dimensions of the irreps ' num2str(group.irrep_dims)]);
        
    else
        
        % Remind the user that symmetry is not operational
        report(spin_system,'no symmetry information available.');
        
    end
    
    % Run the SALC procedure
    if exist('group','var')
        
        % Preallocate the permutation table
        permutation_table=zeros(size(spin_system.bas.basis,1),group.order);
        
        % Compute the permutation table
        parfor n=1:group.order %#ok<*PFBNS>
            group_element=1:spin_system.comp.nspins;
            group_element(spins)=group_element(spins(group.elements(n,:)));
            permuted_basis=spin_system.bas.basis(:,group_element);
            [~,index]=sortrows(permuted_basis);
            permutation_table(:,n)=index;
        end
        
        % Compute irreducible representations
        if spin_system.comp.sym_a1g_only
            
            % Inform the user
            report(spin_system,'Liouville space symmetry mode - fully symmetric irrep only.');
            
            % Prune the permutation table
            symmetry_related_states=unique(sort(permutation_table,2,'ascend'),'rows');
            dimension=size(symmetry_related_states,1);
            
            % Populate the coefficient matrix
            index=unique([kron(ones(group.order,1),(1:dimension)') symmetry_related_states(:) ones(dimension*group.order,1)],'rows');
            coeff_matrix=sparse(index(:,1),index(:,2),index(:,3),dimension,size(spin_system.bas.basis,1));
            
            % Normalize the coefficient matrix
            norms=sqrt(sum(coeff_matrix.^2,2));
            coeff_matrix=spdiags(norms.^(-1),0,dimension,dimension)*coeff_matrix;
            
            % Report back to the user
            report(spin_system,['A1g irrep, ' num2str(dimension) ' states.']);
            
            % Return the projector and dimension
            spin_system.bas.irrep.projector=coeff_matrix';
            spin_system.bas.irrep.dimension=dimension;
            
        else
            
            % Inform the user
            report(spin_system,'WARNING - full symmetry treatment, all irreps will be included.');
            
            % Determine the problem dimension
            basis_dim=size(spin_system.bas.basis,1);
            
            % Loop over irreducible representations
            for n=1:group.n_irreps
                
                % Inform the user
                report(spin_system,['processing irrep #' num2str(n) ' out of ' num2str(group.n_irreps) '...']); 
                
                % Preallocate the coefficient matrix
                coeff_matrix=spalloc(basis_dim,basis_dim,basis_dim*group.order);
                
                % Compute the coefficient matrix
                for k=1:basis_dim
                    
                    % Populate the coefficient matrix
                    for m=1:group.order
                        coeff_matrix(k,permutation_table(k,m))=coeff_matrix(k,permutation_table(k,m))+group.characters(n,m); %#ok<SPRIX>
                    end
                    
                    % Unify the signs
                    non_zero_elements=nonzeros(coeff_matrix(k,:));
                    if any(non_zero_elements)
                        coeff_matrix(k,:)=coeff_matrix(k,:)*sign(non_zero_elements(1)); %#ok<SPRIX>
                    end
                    
                end
                
                % Remove identical rows and all-zero rows
                coeff_matrix=unique(coeff_matrix,'rows');
                coeff_matrix(sum(abs(coeff_matrix),2)<spin_system.tols.liouv_zero,:)=[];
                
                % If the irrep is 2+ dimensional, orthonormalize the SALCs. Otherwise, just normalize the SALCs.
                if group.irrep_dims(n)>1
                    report(spin_system,['WARNING: ' num2str(group.irrep_dims(n)) '-dimensional irrep encountered. Will have to orthogonalize (slow).']);
                    coeff_matrix=orth(full(coeff_matrix'))';
                else
                    for k=1:size(coeff_matrix,1)
                        coeff_matrix(k,:)=coeff_matrix(k,:)/norm(coeff_matrix(k,:));
                    end
                end
                
                % Inform the user and write the irrep into the data structure
                report(spin_system,['irreducible representation #' num2str(n) ', ' num2str(size(coeff_matrix,1)) ' states.']);
                spin_system.bas.irrep(n).projector=coeff_matrix';
                spin_system.bas.irrep(n).dimension=size(coeff_matrix,1);
                
            end
            
        end
        
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,bas)

% Check symmetry parameters
if isfield(bas,'sym_group')
    
    % Check the type
    if (~iscell(bas.sym_group))||any(~cellfun(@ischar,bas.sym_group))
        error('bas.sym_group must be a cell array of strings.');
    end
    
    % Check that bas.sym_spins exists
    if ~isfield(bas,'sym_spins')
        error('bas.sym_spins must be specified alongside bas.sym_group.');
    end
    
    % Check the type
    if (~iscell(bas.sym_spins))||any(~cellfun(@isnumeric,bas.sym_spins))
        error('bas.sym_spins must be a cell array of vectors.');
    end
    
    % Check the dimensions
    if numel(bas.sym_spins)~=numel(bas.sym_group)
        error('bas.sym_group and bas.sym_spins arrays must have the same number of elements.');
    end
    
    % Check the spin indices
    for m=1:length(bas.sym_spins)
        
        % Check for sense
        if any(bas.sym_spins{m}>spin_system.comp.nspins)||any(bas.sym_spins{m}<1)||(numel(bas.sym_spins{m})<2)
            error('incorrect spin labels in bas.sym_spins.');
        end
        
        % Check for intersections
        for n=1:length(bas.sym_spins)
            if (n~=m)&&(~isempty(intersect(bas.sym_spins{m},bas.sym_spins{n})))
                error('same spin is listed in multiple symmetry groups in bas.sym_spins.');
            end
        end
        
    end
    
    % Check the group names
    for n=1:length(bas.sym_group)
        if ~ismember(bas.sym_group{n},{'S2','S3','S4','S4A','S5','S6'})
            error('the group requested in bas.sym_group is not available.');
        end
    end
    
    % Check the irrep switch
    if isfield(bas,'sym_a1g_only')&&(~isnumeric(bas.sym_a1g_only))&&(~islogical(bas.sym_a1g_only))&&((bas.sym_a1g_only~=1)||(bas.sym_a1g_only~=0))
        error('the allowed values for bas.sym_a1g_only are 0 and 1.');
    end
    
else
    
    % Enforce no sys.sym_spins without sys.sym_group
    if isfield(bas,'sym_spins')
        error('bas.sym_group must be specified alongside bas.sym_spins.');
    end
    
    % Enforce no sys.sym_a1g_only without sys.sym_group
    if isfield(bas,'sym_a1g_only')
        error('bas.sym_group must be specified alongside bas.sym_a1g_only.');
    end
    
end

end

% I am regularly asked what an average Internet user can
% do to ensure his security. My first answer is usually
% "nothing, you're screwed".
%
% Bruce Schneier

