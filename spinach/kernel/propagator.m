% Calculates exponential propagators. Syntax:
%
%            P=propagator(spin_system,L,timestep)
%
% returns exp(-i*L*t). The following calculation methods are
% supported:
%
%   'cpu'       - Taylor series with scaling and squaring
%                 on CPU, spmd parallel if possible
%
%   'gpu'       - Taylor series with scaling and squaring
%                 on GPU
%
% The propagator calculation method is chosen by setting the
% sys.enable parameter during the call to create.m function.
%
% Notes: we did have Chebyshev and Newton series here at one
%        point, as well as the Pade method. None of them has 
%        lived up to their marketing.
%
% i.kuprov@soton.ac.uk

function P=propagator(spin_system,L,timestep)

% Check consistency
grumble(L,timestep);

% Inform the user about matrix densities
report(spin_system,['Liouvillian dimension ' num2str(size(L,1)) ', nnz ' num2str(nnz(L)) ...
                    ', density ' num2str(100*nnz(L)/numel(L)) '%, sparsity ' num2str(issparse(L))]);

% Set a shorthand for -i*L*dt
A=-1i*L*timestep;

% Fast bypass for small matrices
if size(A,1)<spin_system.tols.small_matrix
    P=expm(full(A)); return;
end

% Check the cache
if ismember('caching',spin_system.sys.enable)&&(size(A,1)>500)
    
    % Generate the cache record name
    filename=[spin_system.sys.scratch '/p_' md5_hash(A) '.mat'];
    
    % Try loading the cache record
    if exist(filename,'file')
        report(spin_system,'cache record found and used.'); load(filename,'P'); 
        report(spin_system,['propagator dimension ' num2str(size(P,1)) ', nnz ' num2str(nnz(P))...
                            ', density ' num2str(100*nnz(P)/numel(P)) '%, sparsity ' num2str(issparse(P))]); %#ok<NODEF>
        return;
    else
        report(spin_system,'cache record not found, computing from scratch...');
    end
    
end

% Get the norm
mat_norm=norm(A,1);

% Warn the user if the norm is too big
if mat_norm>1024
    
    % If the user is really pushing it, take precautionary measures
    report(spin_system,'WARNING - the time step requested greatly exceeds the timescale of system dynamics.');
    report(spin_system,'WARNING - exponentiation tolerance will be set to 1e-14.');
    spin_system.tols.prop_chop=1e-14;
    
elseif mat_norm>16
    
    % Inform the user just in case
    report(spin_system,'WARNING - the time step requested exceeds the timescale of system dynamics.');
    
end

% Determine scaling and squaring parameters
n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings;
report(spin_system,['scaling by ' num2str(scaling_factor) ' and squaring ' num2str(n_squarings) ' times.']);

% Scale and clean up the matrix
if scaling_factor>1, A=A/scaling_factor; end
A=clean_up(spin_system,A,spin_system.tols.prop_chop);

% Get the propagator
if ismember('gpu',spin_system.sys.enable)&&(size(A,1)>500)
    
    % Run Taylor series procedure on the GPU
    A=gpuArray(A); P=speye(size(A));
    next_term=gpuArray.speye(size(A)); n=1;
    while nnz(next_term)>0
        
        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end
        
        % Eliminate small elements
        next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
        
        % Add to the total and increment the counter
        P=P+gather(next_term); n=n+1;
        
    end
    
    % Inform the user
    report(spin_system,['Taylor series converged on GPU in ' num2str(n) ' iterations.']);
    
else
    
    % Run Taylor series procedure on the CPU
    P=speye(size(A)); next_term=P; n=1;
    while nnz(next_term)>0
        
        % Compute the next term
        if issparse(A)
            next_term=(1/n)*A*next_term;
        else
            next_term=(1/n)*next_term*A;
        end
        
        % Eliminate small elements
        next_term=clean_up(spin_system,next_term,spin_system.tols.prop_chop);
        
        % Add to the total and increment the counter
        P=P+next_term; n=n+1;
        
    end
    
    % Inform the user
    report(spin_system,['Taylor series converged on CPU in ' num2str(n) ' iterations.']);
    
end

% Reclaim memory
clear('A','next_term');

% Clean up the result
P=clean_up(spin_system,P,spin_system.tols.prop_chop);

% Run the squaring stage
if n_squarings>0
    
    % Run the appropriate squaring process
    if ismember('gpu',spin_system.sys.enable)&&(size(P,1)>500)
        
        % Move the array to GPU
        P=gpuArray(P);
        
        % Run the squaring on the GPU
        for n=1:n_squarings
            
            % Inform the user
            report(spin_system,['GPU squaring step ' num2str(n) '...']);
            
            % Square the propagator
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
            
        end
        
    elseif (~isworkernode)&&(nnz(P)>1e6)&&issparse(P)
        
        % Run codistributed CPU squaring
        spmd
            
            % Codistribute the propagator
            P=codistributed(P,codistributor2dbc());
            
            % Run the squaring stage
            for n=1:n_squarings
                
                % Inform the user
                report(spin_system,['codistributed CPU squaring step ' num2str(n) '...']);
                
                % Square the propagator
                P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
                
            end
            
        end
        
    else
        
        % Run serial CPU squaring
        for n=1:n_squarings
            
            % Inform the user
            report(spin_system,['CPU squaring step ' num2str(n) '...']);
            
            % Square the propagator
            P=clean_up(spin_system,P*P,spin_system.tols.prop_chop);
            
        end
      
    end
    
end

% Gather the propagator
P=gather(P);

% Inform the user
report(spin_system,['propagator dimension ' num2str(size(P,1)) ', nnz ' num2str(nnz(P))...
                    ', density ' num2str(100*nnz(P)/numel(P)) '%, sparsity ' num2str(issparse(P))]);

% Write the cache record
if ismember('caching',spin_system.sys.enable)&&any(size(P)>500)
    
    % Save the propagator
    save(filename,'P'); report(spin_system,'cache record saved.');
    
end

end

% Consistency enforcement
function grumble(L,timestep)
if (~isnumeric(L))||(~isnumeric(timestep))
    error('both L and timestep must be numeric.');
end
if size(L,1)~=size(L,2)
    error('L argument must be a square matrix.');
end
if ~isscalar(timestep)
    error('timestep must be a scalar.');
end
end

% To preserve one's mind intact through a modern college education is a
% test of courage and endurance, but the battle is worth it and the sta-
% kes are the highest possible to man: the survival of reason.
%
% Ayn Rand

