% Krylov propagation function. Avoids matrix exponentiation, but can be 
% slow. Should be used when the Liouvillian exponential does not fit in-
% to the system memory, but the Liouvillian itself does. Syntax:
%
%     answer=krylov(spin_system,L,coil,rho,time_step,nsteps,output)
%
% Arguments for Liouville space calculations:
%
%      L      - the Liouvillian to be used during evolution
%
%      rho    - the initial state vector or a horizontal stack thereof
%
%      output - a string giving the type of evolution that is required
%
%                'final' - returns the final state vector or a horizontal
%                          stack thereof.
%
%                'trajectory' - returns the stack of state vectors giving
%                               the trajectory of the system starting from
%                               rho with the user-specified number of steps
%                               and step length.
%
%                'total'   - returns the integral of the observable trace
%                            from the simulation start to infinity. This
%                            option requires the presence of relaxation.
%
%                'refocus' - evolves the first vector for zero steps,
%                            second vector for one step, third vector for
%                            two steps, etc., consistent with the second
%                            stage of evolution in the indirect dimension
%                            after a refocusing pulse.
%
%                'observable' - returns the time dynamics of an observable
%                               as a vector (if starting from a single ini-
%                               tial state) or a matrix (if starting from a
%                               stack of initial states).
%
%                'multichannel' - returns the time dynamics of several
%                                 observables as rows of a matrix. Note
%                                 that destination state screening may be
%                                 less efficient when there are multiple
%                                 destinations to screen against.
%
%      coil   - the detection state, used when 'observable' is specified as
%               the output option. If 'multichannel' is selected, the coil
%               should contain multiple columns corresponding to individual
%               observable vectors.
%
% This function does not support Hilber space formalisms.
%
% i.kuprov@soton.ac.uk

function answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output)

% Check consistency
grumble(L,coil,rho,timestep,nsteps,output);

% Inform the user about matrix densities
report(spin_system,['Liouvillian dimension ' num2str(size(L,1)) ', nnz ' num2str(nnz(L)) ...
                    ', density ' num2str(100*nnz(L)/numel(L)) '%, sparsity ' num2str(issparse(L))]);

% Estimate the norm of i*L*dt
norm_mat=norm(L,1)*abs(timestep);

% Inform the user about the norm
report(spin_system,['norm of i*L*dt: ' num2str(norm_mat)]);

% Determine the number of substeps
nsubsteps=ceil(norm_mat);

% Inform the user about the schedule
report(spin_system,['taking ' num2str(nsteps) ' Krylov steps with '...
                              num2str(nsubsteps) ' substeps each.']);

% Warn about silly calls
if nsubsteps>100
    report(spin_system,'WARNING - the number of substeps is large, consider using evolution() here.');
end

% Estimate the scaling coefficient
scaling=max([1 norm(rho,1)]);

% Scale the vector
rho=rho/scaling;

% Compute the generator matrix
A=-1i*L*(timestep/nsubsteps);

% Upload data to GPU or optimise layout
if ismember('gpu',spin_system.sys.enable)
    A=gpuArray(A); rho=gpuArray(full(rho));
    coil=gpuArray(coil); location='GPU';
else
    location='CPU'; rho=full(rho);
end

% Decide the output type
switch output
    
    case 'final'
        
        % Loop over steps
        for n=1:nsteps
            
            % Inform the user
            report(spin_system,[location ' Krylov step ' num2str(n) ' out of ' num2str(nsteps) '...']);
            
            % Loop over substeps
            for k=1:nsubsteps
                
                % Taylor series
                next_term=rho; m=1;
                while nnz(abs(next_term)>eps)>0
                    next_term=(A*next_term)*(1/m);
                    rho=rho+next_term; m=m+1;
                end
                
            end
            
        end
        
        % Assign the answer
        answer=gather(rho);
                
    case 'trajectory'
        
        % Preallocate the answer and set the starting point
        answer=zeros(size(rho,1),nsteps+1);
        answer(:,1)=gather(rho);
        
        % Loop over steps
        for n=1:nsteps
            
            % Inform the user
            report(spin_system,[location ' Krylov step ' num2str(n) ' out of ' num2str(nsteps) '...']);
            
            % Loop over substeps
            for k=1:nsubsteps
                
                % Taylor series
                next_term=rho; m=1;
                while nnz(abs(next_term)>eps)>0
                    next_term=(A*next_term)*(1/m);
                    rho=rho+next_term; m=m+1;
                end
                
            end
            
            % Assign the answer
            answer(:,n+1)=gather(rho);
            
        end
        
    case 'refocus'
        
        % Loop over steps
        for n=2:size(rho,2)
            
            % Inform the user
            report(spin_system,[location ' Krylov step ' num2str(n) ' out of ' num2str(nsteps) '...']);
            
            % Loop over substeps
            for k=1:nsubsteps
                
                % Taylor series
                next_term=rho(:,n:end); m=1;
                while nnz(abs(next_term)>eps)>0
                    next_term=(A*next_term)*(1/m); 
                    rho(:,n:end)=rho(:,n:end)+next_term; m=m+1;
                end
                
            end
            
        end
        
        % Assign the answer
        answer=gather(rho);
        
    case 'observable'
        
        % Preallocate the answer
        answer=zeros(nsteps+1,size(rho,2));
        
        % Set the initial point
        answer(1,:)=gather(coil'*rho);
        
        % Loop over steps
        for n=1:nsteps
            
            % Inform the user
            report(spin_system,[location ' Krylov step ' num2str(n) ' out of ' num2str(nsteps) '...']);
             
            % Loop over substeps
            for k=1:nsubsteps
                
                % Taylor series
                next_term=rho; m=1;
                while nnz(abs(next_term)>eps)>0
                    next_term=(A*next_term)*(1/m);
                    rho=rho+next_term; m=m+1;
                end
                
            end
            
            % Assign the answer
            answer(n+1,:)=gather(coil'*rho);
            
        end
        
    case 'multichannel'
        
        % Preallocate the answer
        answer=zeros(size(coil,2),nsteps+1);
        
        % Set the initial point
        answer(:,1)=answer(:,1)+gather(coil'*rho);
        
        % Loop over steps
        for n=1:nsteps
            
            % Inform the user
            report(spin_system,[location ' Krylov step ' num2str(n) ' out of ' num2str(nsteps) '...']);
            
            % Loop over substeps
            for k=1:nsubsteps
                
                % Taylor series
                next_term=rho; m=1;
                while nnz(abs(next_term)>eps)>0
                    next_term=(A*next_term)*(1/m);
                    rho=rho+next_term; m=m+1;
                end
                
            end
            
            % Assign the answer
            answer(:,n+1)=gather(coil'*rho);
            
        end
        
    otherwise
        
        error('unknown output option.');
        
end

end

% Consistency enforcement
function grumble(L,coil,rho,timestep,nsteps,output)
if ~isnumeric(L)
    error('Liouvillian must be numeric.');
end
if ~isnumeric(coil)
    error('coil argument must be numeric.');
end
if (~isnumeric(rho))&&(~iscell(rho))
    error('rho argument must either be numeric or a cell array');
end
if ~isnumeric(timestep)
    error('timestep argument must be numeric.');
end
if ~isnumeric(nsteps)
    error('nsteps must be numeric.');
end
if (~ischar(output))||(~ismember(output,{'observable','final',...
   'trajectory','total','multichannel','refocus'}))
    error('observable argument must be a valid character string.');
end
end

% Health and Safety were a pair of rats that lived in
% the shrubbery in front of the Southampton Chemistry
% building. The little beasts have died of old age. 

