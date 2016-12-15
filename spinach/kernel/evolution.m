% Time evolution function. Performs all types of time propagation with
% automatic trajectory level state space restriction. Syntax:
%
%         answer=evolution(spin_system,L,coil,rho,timestep,...
%                          nsteps,output,destination)
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
%      destination - (optional) the state to be used for destination state
%                    screening.
%
% Arguments for Hilbert space calculations:
%
%       L         - Hamiltonian matrix
%
%       coil      - observable operator (if any)
%
%       rho       - initial density matrix
%
%       timestep  - duration of a single time step (seconds)
%
%       nsteps    - number of steps to take
%
%       output    - a string giving the type of evolution that is required
%
%                'final' - returns the final density matrix.
%
%                'trajectory' - returns a cell array of density matrices
%                               giving the trajectory of the system star-
%                               ting from rho with the user-specified num-
%                               ber of steps and step length.
%
%                'refocus' - evolves the first matrix for zero steps,
%                            second matrix for one step, third matrix for
%                            two steps, etc., consistent with the second
%                            stage of evolution in the indirect dimension
%                            after a refocusing pulse.
%
%                'observable' - returns the time dynamics of an observable
%                               as a vector.
%
%       destination - this argument is ignored.
%
% Calculation of final states and observables in Hilbert space is parallel-
% ized and tested all the way to 128-core (16 nodes, 8 cores each) configu-
% rations. Parallelization of the trajectory calculation does not appear to
% yield any benefits due to large amount of inter-thread communication.
%
% See http://dx.doi.org/10.1063/1.3679656 for further information.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk
% ohad.levinkron@weizmann.ac.il

function answer=evolution(spin_system,L,coil,rho,timestep,nsteps,output,destination)

% Check consistency
grumble(L,coil,rho,timestep,nsteps,output);

% Decide how to proceed
switch spin_system.bas.formalism
    
    case {'sphten-liouv','zeeman-liouv'}
        
        % Decide whether to normalize the coils
        if ~isempty(coil)
            if ismember('norm_coil',spin_system.sys.disable)
                report(spin_system,'coils have not been normalized.');
            else
                for n=1:size(coil,2)
                    coil(:,n)=coil(:,n)/norm(coil(:,n));
                end
                report(spin_system,'coils have been normalized.');
            end
        end
        
        % Apply trajectory-level reduction algorithms
        report(spin_system,'trying to reduce the problem dimension...');
        
        % Decide the screening algorithm
        if ismember('dss',spin_system.sys.disable)
            
            % If DSS is disabled, run forward reduction
            report(spin_system,'WARNING - destination state screening is disabled.');
            projectors=reduce(spin_system,L,rho);
            
        elseif exist('destination','var')&&(~isempty(destination))
            
            % If destination state is supplied, use it for DSS
            report(spin_system,'destination state screening using supplied state.');
            projectors=reduce(spin_system,L',destination);
            
        elseif ismember(output,{'observable'})
            
            % If coil state is supplied, use it for DSS
            report(spin_system,'destination state screening using coil state.');
            projectors=reduce(spin_system,L',coil);
            
        else
            
            % Default to the usual forward screening
            report(spin_system,'destination state screening is not applicable.');
            projectors=reduce(spin_system,L,rho);
            
        end
        
        % Run the evolution
        switch output
            
            case 'final'
                
                % Create arrays of projections
                nsubs=numel(projectors); L_sub=cell(nsubs,1);
                rho_sub=cell(nsubs,1); local_answers=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                end
                
                % Loop in parallel over independent subspaces
                parfor sub=1:nsubs
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub};
                    
                    % Propagate the system
                    if ((nnz(L_loc)<spin_system.tols.krylov_switchover)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Move to GPU or optimise layout
                        if ismember('gpu',spin_system.sys.enable)
                            rho_loc=gpuArray(full(rho_loc)); P=gpuArray(P);
                            report(spin_system,'propagating the system on GPU...');
                        else
                            rho_loc=full(rho_loc);
                            report(spin_system,'propagating the system on CPU...');
                        end
                        
                        % Propagate the system
                        for n=1:nsteps
                            rho_loc=P*rho_loc;
                        end
                        
                        % Gather the answer
                        rho_loc=gather(rho_loc);
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        rho_loc=krylov(spin_system,L_loc,[],rho_loc,timestep,nsteps,'final');
                        
                    end
                    
                    % Project back into the full space
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    local_answers{sub}=proj_loc*sparse(rho_loc);
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; proj_loc=[]; P=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
                % Sum up the results (no clean-up needed)
                report(spin_system,'gathering data from the nodes...');
                answer=spalloc(size(rho,1),size(rho,2),0);
                for sub=1:numel(projectors)
                    answer=answer+local_answers{sub};
                    local_answers{sub}=[];
                end
                report(spin_system,'data retrieval finished.');
                
            case 'total'
                
                % Create arrays of projections
                nsubs=numel(projectors); L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=full(projectors{sub}'*rho);
                    coil_sub{sub}=projectors{sub}'*coil;
                end
                
                % Start with zero
                answer=0;
                
                % Loop in parallel over independent subspaces
                parfor sub=1:numel(projectors)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Compute the integral
                    report(spin_system,'integrating the trajectory...');
                    answer_loc=real(coil_sub{sub}'*((1i*L_sub{sub})\rho_sub{sub}));
                    
                    % Add the subspace to the total
                    answer=answer+answer_loc;
                    
                    % Update the user
                    report(spin_system,'integration finished.');
                    
                end
                
                % Make sure a full array is returned
                answer=full(answer);
                
            case 'trajectory'
                
                % Create arrays of projections
                nsubs=numel(projectors); L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); local_answers=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                end
                
                % Loop in parallel over independent subspaces
                parfor sub=1:nsubs
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub};
                    
                    % Propagate the system
                    if ((nnz(L_loc)<spin_system.tols.krylov_switchover)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                    
                        % Move to GPU or optimise layout
                        if ismember('gpu',spin_system.sys.enable)
                            rho_loc=gpuArray(full(rho_loc)); P=gpuArray(P);
                            answer=zeros(size(rho_loc,1),nsteps+1,'gpuArray');
                            report(spin_system,'propagating the system on GPU...');
                        else
                            rho_loc=full(rho_loc);
                            answer=zeros(size(rho_loc,1),nsteps+1);
                            report(spin_system,'propagating the system on CPU...');
                        end

                        % Set the starting point
                        answer(:,1)=rho_loc;
                        
                        % Propagate the system
                        for n=1:nsteps
                            answer(:,n+1)=P*answer(:,n);
                        end
                        
                        % Gather the answer
                        answer=gather(answer);
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer=krylov(spin_system,L_loc,[],rho_loc,timestep,nsteps,'trajectory')
                        
                    end
                    
                    % Project back into the full space
                    answer=clean_up(spin_system,answer,spin_system.tols.zte_tol);
                    local_answers{sub}=proj_loc*sparse(answer);
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; proj_loc=[]; P=[]; answer=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
                % Sum up the results (no clean-up needed)
                report(spin_system,'gathering data from the nodes...');
                answer=spalloc(size(rho,1),nsteps+1,0);
                for sub=1:numel(projectors)
                    answer=answer+local_answers{sub};
                    local_answers{sub}=[];
                end
                report(spin_system,'data retrieval finished.');
                
            case 'refocus'
                
                % Create arrays of projections
                nsubs=numel(projectors); L_sub=cell(nsubs,1);
                rho_sub=cell(nsubs,1); local_answers=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                end
                
                % Loop over independent subspaces
                parfor sub=1:nsubs
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; proj_loc=projectors{sub};
                    
                    % Propagate the system
                    if ((nnz(L_loc)<spin_system.tols.krylov_switchover)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Move to GPU or optimise layout
                        if ismember('gpu',spin_system.sys.enable)
                            rho_loc=gpuArray(full(rho_loc)); P=gpuArray(P);
                            report(spin_system,'propagating the system on GPU...');
                        else
                            rho_loc=full(rho_loc);
                            report(spin_system,'propagating the system on CPU...');
                        end
                        
                        % Propagate the system
                        for n=2:size(rho_loc,2)
                            rho_loc(:,n:end)=P*rho_loc(:,n:end);
                        end
                        
                        % Gather the answer
                        rho_loc=gather(rho_loc);
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        rho_loc=krylov(spin_system,L_loc,[],rho_loc,timestep,'refocus');
                        
                    end
                    
                    % Project back into the full space
                    rho_loc=clean_up(spin_system,rho_loc,spin_system.tols.zte_tol);
                    local_answers{sub}=proj_loc*sparse(rho_loc);
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; proj_loc=[]; P=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
                % Sum up the results (no clean-up needed)
                report(spin_system,'gathering data from the nodes...');
                answer=spalloc(size(rho,1),nsteps+1,0);
                for sub=1:numel(projectors)
                    answer=answer+local_answers{sub};
                    local_answers{sub}=[];
                end
                report(spin_system,'data retrieval finished.');
                
            case 'observable'
                
                % Create arrays of projections
                nsubs=numel(projectors); L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                    coil_sub{sub}=projectors{sub}'*coil;
                end
                
                % Preallocate the answer
                answer=zeros(nsteps+1,size(rho,2));
                
                % Loop over independent subspaces
                parfor sub=1:nsubs
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(nsubs) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; coil_loc=coil_sub{sub};
                    
                    % Preallocate the local answer
                    answer_loc=zeros(nsteps+1,size(rho_loc,2));
                    
                    % The first point does not require propagation
                    answer_loc(1,:)=answer_loc(1,:)+coil_loc'*rho_loc;
                    
                    % Propagate the system
                    if ((nnz(L_loc)<spin_system.tols.krylov_switchover)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Move to GPU or optimise layout
                        if ismember('gpu',spin_system.sys.enable)
                            rho_loc=gpuArray(full(rho_loc));
                            coil_loc=gpuArray(full(coil_loc)); P=gpuArray(P);
                            report(spin_system,'propagating the system on GPU...');
                        else
                            rho_loc=full(rho_loc); coil_loc=full(coil_loc);
                            report(spin_system,'propagating the system on CPU...');
                        end
                        
                        % Propagate the system
                        for n=1:nsteps
                            rho_loc=P*rho_loc;
                            answer_loc(n+1,:)=answer_loc(n+1,:)+gather(coil_loc'*rho_loc);
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_loc=krylov(spin_system,L_loc,coil_loc,rho_loc,timestep,nsteps,'observable');
                        
                    end
                    
                    % Add to the total
                    answer=answer+answer_loc;
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; P=[]; answer_loc=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
                % Make sure a full array is returned
                answer=full(answer);
                
            case 'multichannel'
                
                % Create arrays of projections
                nsubs=numel(projectors); L_sub=cell(nsubs,1); 
                rho_sub=cell(nsubs,1); coil_sub=cell(nsubs,1);
                report(spin_system,'splitting the space...');
                for sub=1:nsubs
                    L_sub{sub}=projectors{sub}'*L*projectors{sub};
                    rho_sub{sub}=projectors{sub}'*rho;
                    coil_sub{sub}=projectors{sub}'*coil;
                end
                
                % Preallocate the answer
                answer=zeros(size(coil,2),nsteps+1);
                
                % Loop over independent subspaces
                parfor sub=1:length(projectors)
                    
                    % Inform the user
                    report(spin_system,['evolving subspace ' num2str(sub) ' of ' num2str(numel(projectors)) '...']);
                    
                    % Grab local copies
                    L_loc=L_sub{sub}; rho_loc=rho_sub{sub}; coil_loc=coil_sub{sub};
                    
                    % Propagate the system
                    if ((nnz(L_loc)<spin_system.tols.krylov_switchover)||...
                        ismember('krylov',spin_system.sys.disable))&&...
                      (~ismember('krylov',spin_system.sys.enable))
                        
                        % Get the exponential propagator
                        report(spin_system,'building the propagator...');
                        P=propagator(spin_system,L_loc,timestep);
                        
                        % Preallocate the local answer
                        answer_loc=zeros(size(coil_loc,2),nsteps+1);
                        
                        % The first point does not require propagation
                        answer_loc(:,1)=answer_loc(:,1)+coil_loc'*rho_loc;
                        
                        % Move to GPU or optimise layout
                        if ismember('gpu',spin_system.sys.enable)
                            rho_loc=gpuArray(full(rho_loc));
                            coil_loc=gpuArray(full(coil_loc)); P=gpuArray(P);
                            report(spin_system,'propagating the system on GPU...');
                        else
                            rho_loc=full(rho_loc); coil_loc=full(coil_loc);
                            report(spin_system,'propagating the system on CPU...');
                        end
                        
                        % Propagate the system
                        for n=1:nsteps
                            rho_loc=P*rho_loc;
                            answer_loc(:,n+1)=answer_loc(:,n+1)+gather(coil_loc'*rho_loc);
                        end
                        
                    else
                        
                        % For very large subspaces use Krylov propagation
                        report(spin_system,'large Liouvillian, propagating using Krylov algorithm... ');
                        answer_loc=krylov(spin_system,L_loc,coil_loc,rho_loc,timestep,nsteps,'multichannel');
                        
                    end
                    
                    % Add to the total
                    answer=answer+answer_loc;
                    
                    % Deallocate memory in a transparent way
                    rho_loc=[]; L_loc=[]; P=[]; answer_loc=[]; %#ok<NASGU>
                    
                    % Update the user
                    report(spin_system,'propagation finished.');
                    
                end
                
                % Make sure a full array is returned
                answer=full(answer);
                
            otherwise
                
                error('invalid output option.');
                
        end
        
    case 'zeeman-hilb'
        
        % Decide whether to normalize the coils
        if ~isempty(coil)
            if ismember('norm_coil',spin_system.sys.disable)
                report(spin_system,'coil has not been normalized.');
            else
                coil=coil/norm(coil,'fro');
                report(spin_system,'coil has been normalized.');
            end
        end
        
        % Run the evolution
        switch output
            
            case 'final'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,timestep);
                        
                        % Report to the user
                        report(spin_system,'propagating the system...');
                        
                        % Decide parallelization style
                        if iscell(rho)
                            
                            % Run in parallel over cells
                            parfor k=1:numel(rho)
                                for n=1:nsteps
                                    rho{k}=P*rho{k}*P'
                                end
                            end
                            
                            % Return the result
                            answer=rho;
                            
                        else
                            
                            % If inside a parallel loop, avoid spmd
                            if isworkernode
                                
                                % Propagate vectors and covectors
                                for n=1:nsteps
                                    rho=P*rho*P';
                                end
                                
                                % Return the result
                                answer=rho;
                                
                            else
                                
                                % Create covector array
                                cov=speye(size(rho));
                                
                                % Run in parallel over columns
                                spmd
                                    
                                    % Slice density matrix and covectors
                                    codist=codistributor('1d',2);
                                    rho=codistributed(rho,codist);
                                    cov=codistributed(cov,codist);
                                    
                                    % Propagate vectors and covectors
                                    for n=1:nsteps
                                        rho=P*rho; cov=P*cov;
                                    end
                                    
                                end
                                
                                % Return the result
                                answer=gather(rho)*gather(cov)';
                                
                            end
                            
                        end
                        
                        % Report to the user
                        report(spin_system,'propagation finished.');
                        
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
            case 'observable'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,timestep);
                        
                        % Report to the user
                        report(spin_system,'propagating the system...');
                        
                        % If a stack is received, run 2D acquisition
                        if iscell(rho)
                            
                            % Preallocate the answer
                            answer=zeros(nsteps+1,numel(rho));
                            
                            % Loop over the elements of the stack
                            parfor k=1:numel(rho)
                                
                                % Grab a local copy
                                rho_local=rho{k};
                                
                                % Loop over time steps
                                for n=1:(nsteps+1)
                                
                                    % Compute the observable
                                    answer(n,k)=hdot(coil,rho_local);
                                
                                    % Step forward
                                    rho_local=P*rho_local*P';
                                    
                                end
                                
                            end
                        
                        % If inside a parallel loop, avoid spmd
                        elseif isworkernode
                            
                            % Preallocate the answer
                            answer=zeros(nsteps+1,1);                   
                            
                            % Loop over time steps
                            for n=1:(nsteps+1)
                                
                                % Compute the observable
                                answer(n)=hdot(coil,rho);
                                
                                % Step forward
                                rho=P*rho*P';
                                
                            end
                            
                        % Otherwise use Kuprov-Edwards parallel split
                        else
                            
                            % Create covector array
                            cov=speye(size(rho));
                            
                            % Parallel processing
                            spmd
                                
                                % Slice density matrix and covectors
                                codist=codistributor('1d',2);
                                rho=codistributed(rho,codist);
                                cov=codistributed(cov,codist);
                                
                                % Localize the problem
                                rho_local=getLocalPart(rho);
                                cov_local=getLocalPart(cov);
                                fid_local=zeros(nsteps+1,1);
                                
                                % Loop over time steps
                                for n=1:(nsteps+1)
                                    
                                    % Write local fid
                                    fid_local(n)=hdot(coil*cov_local,rho_local);
                                    
                                    % Step forward on rho and cov
                                    rho_local=P*rho_local;
                                    cov_local=P*cov_local;
                                    
                                end
                                
                                % Collect the results
                                answer=gplus(fid_local,1);
                                
                            end
                            
                            % Return the result
                            answer=answer{1};
                        
                        end
                        
                        % Report to the user
                        report(spin_system,'propagation finished.');
                                                     
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
                % Make sure a full array is returned
                answer=full(answer);
                
            case 'trajectory'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,timestep);
                        
                        % Preallocate the answer
                        answer=cell(1,nsteps+1);
                        
                        % Compute the trajectory
                        report(spin_system,'propagating the system...');
                        for n=1:(nsteps+1)
                            answer{n}=rho; rho=P*rho*P';
                        end
                        report(spin_system,'propagation finished.');
                        
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
            case 'refocus'
                
                % Decide how to proceed
                switch spin_system.bas.approximation
                    
                    case 'none'
                        
                        % Compute the exponential propagator
                        P=propagator(spin_system,L,timestep);
                        
                        % Propagate the system
                        report(spin_system,'propagating the system...');
                        parfor k=2:numel(rho)
                            
                            % Grab a local copy
                            rho_local=rho{k};
                            
                            % Propagate local copy
                            for n=2:k
                                rho_local=P*rho_local*P';
                            end
                            
                            % Return local copy
                            rho{k}=rho_local;
                            
                        end
                        answer=rho; report(spin_system,'propagation finished.');
                        
                    otherwise
                        
                        % Complain and bomb out
                        error('unrecognised approximation specification.');
                        
                end
                
            otherwise
                
                % Complain and bomb out
                error('unknown evolution option for zeeman-hilb formalism.');
                
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
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

% Degrees of ability vary, but the basic principle remains the same: the
% degree of a man's independence, initiative and personal love for his work
% determines his talent as a worker and his worth as a man. Independence is
% the only gauge of human virtue and value. What a man is and makes of him-
% self; not what he has or hasn't done for others. There is no substitute
% for personal dignity. There is no standard of personal dignity except
% independence.
%
% Ayn Rand, "The Fountainhead"

