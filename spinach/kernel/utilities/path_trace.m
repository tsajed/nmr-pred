% Liouvillian path tracing. Treats the user-supplied Liouvillian 
% as the adjacency matrix of a graph, computes the weakly connect-
% ed subgraphs of that graph and returns a cell array of project-
% ors into the corresponding independently evolving subspaces.
% Syntax:
%              projectors=reduce(spin_system,L,rho)
%
% where L is the Liouvillian and rho is the initial state (in the
% case of source state screening) or the detection state (if des-
% tination state screening is used). The output is a cell array of
% projectors into independently evolving reduced subspaces. Those
% projectors are to be used as follows:
%
%           L_reduced=P'*L*P;     rho_reduced=P'*rho;
%
% Further information is available here:
%
%         http://link.aip.org/link/doi/10.1063/1.3398146
%         http://dx.doi.org/10.1016/j.jmr.2011.03.010
%
% i.kuprov@soton.ac.uk

function projectors=path_trace(spin_system,L,rho)

% Check the input
grumble(spin_system,L,rho);

% Check run conditions
if ismember('pt',spin_system.sys.disable)
    
    % Return a unit projector if path tracing is disabled
    report(spin_system,'WARNING - path tracing disabled by the user.');
    projectors={speye(size(L))}; return
    
elseif size(L,2)<256
    
    % Return a unit projector if the space is small anyway
    report(spin_system,'small space - path tracing skipped.');
    projectors={speye(size(L))}; return
    
else
    
    % Report to the user
    report(spin_system,['analyzing ' num2str(size(L,1)) '-dimensional state space.']);
    report(spin_system,['Liouvillian zero tolerance ' num2str(spin_system.tols.liouv_zero)]);
    report(spin_system,['cross-term tolerance for path drop ' num2str(spin_system.tols.path_drop)]);
    report(spin_system,['population tolerance for subspace drop ' num2str(spin_system.tols.subs_drop)]);
    
end

% Get the connectivity matrix
G=(abs(L)>spin_system.tols.path_drop);

% Make sure isolated states do not get lost
G=or(G,transpose(G)); G=or(G,speye(size(G)));

% Get the weakly connected subgraphs
member_states=scomponents(G);

% Determine the number of subspaces
n_subspaces=max(member_states);
report(spin_system,['found ' num2str(n_subspaces) ' non-interacting subspaces.']);
report(spin_system,'running subspace population analysis...');

% Analyze independent subspaces
subspace_important=true(n_subspaces,1);
tolerance=spin_system.tols.subs_drop;
parfor n=1:n_subspaces
    
    % Determine the importance
    subspace_important(n)=(norm(rho.*(member_states==n),1)>tolerance);
    
end
significant_subspaces=find(subspace_important);
n_subspaces=numel(significant_subspaces);

% Preallocate projectors and counters
projectors=cell(1,n_subspaces);

% Build projectors into significant subspaces
parfor n=1:n_subspaces
    
    % Find states populating the current subspace
    state_index=find(member_states==significant_subspaces(n));

    % Determine subspace dimension
    subspace_dim=numel(state_index);
    
    % Build the projector into the current subspace
    projectors{n}=sparse(state_index,1:subspace_dim,ones(1,subspace_dim),size(L,1),subspace_dim);
                     
end

% Report to the user
for n=1:n_subspaces
    report(spin_system,['populated subspace found, dimension='  num2str(size(projectors{n},2))]);
end
report(spin_system,['keeping a total of ' num2str(n_subspaces) ' independent subspace(s) '...
                    'of total dimension ' num2str(sum(cellfun(@(x)size(x,2),projectors)))]);
                
% Merge small subspaces
if ismember('merge',spin_system.sys.disable)
    
    % Inform the user
    report(spin_system,'WARNING - small subspace merging disabled by the user.');

else
    
    % Inform the user
    report(spin_system,'merging small subspaces...');
    
    % Compile dimension statistics
    subspace_dims=cellfun(@(x)size(x,2),projectors);
    
    % Call the bin packer
    bins=binpack(subspace_dims,spin_system.tols.merge_dim);
    
    % Group the subspaces
    new_projectors=cell(numel(bins),1);
    for n=1:numel(bins)
        new_projectors{n}=[projectors{bins{n}}];
    end
    projectors=new_projectors;
    
    % Report to the user
    for n=1:numel(projectors)
        report(spin_system,['working subspace ' num2str(n) ', dimension '  num2str(size(projectors{n},2))]);
    end
    
end

end

% Consistency enforcement
function grumble(spin_system,L,rho)
if ~ismember(spin_system.bas.formalism,{'zeeman-liouv','sphten-liouv'})
    error('path tracing is only available for zeeman-liouv and sphten-liouv formalisms.');
end
if (~isnumeric(L))||(~isnumeric(rho))
    error('both inputs must be numeric.');
end
if size(L,1)~=size(L,2)
    error('Liouvillian must be square.');
end
if size(L,2)~=size(rho,1)
    error('Liouvillian and state vector dimensions must be consistent.');
end
end

% "My dear fellow, who will let you?"
% "That's not the point. The point is, who will stop me?"
%
% Ayn Rand, "The Fountainhead"

