% Propagation step function. Uses Krylov propagation and sparse exponenti-
% ation where appropriate. Syntax:
%
%                   rho=step(spin_system,L,rho,time_step)
%
% Arguments:
%
%      L          -  Liouvillian or Hamiltonian to be used for propagation
%
%      rho        -  state vector or density matrix to be propagated
%
%      time_step  -  length of the time step to take
%
% Note: the peculiar sequence of algebraic operations in the Krylov proce-
%       dure is designed to minimize the memory footprint.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function rho=step(spin_system,L,rho,time_step)

% Check consistency
grumble(L,time_step);

% Decide how to proceed
switch spin_system.bas.formalism
    
    case {'sphten-liouv','zeeman-liouv'}
        
        % Decide how to proceed
        if ismember('expv',spin_system.sys.disable)||...
           (size(L,1)<spin_system.tols.small_matrix)
            
            % Use Matlab's expm if the user insists
            rho=expm(-1i*L*time_step)*rho;
            
        else
            
            % Estimate the norm of the i*L*dt matrix
            norm_mat=norm(L,1)*abs(time_step);
            
            % Determine the number of time steps
            nsteps=ceil(norm_mat/5);
            
            % Warn the user if the number is too large
            if nsteps>100
                report(spin_system,['WARNING: ' num2str(nsteps) ' substeps required, consider using evolution() here.']);
            end
            
            % Run codistributed if appropriate
            if (size(rho,2)>1)&&(~isworkernode)&&...
               (~ismember('gpu',spin_system.sys.enable))
            
                spmd
                
                    % Create codistributor object
                    codist=codistributor1d(2);
                    
                    % Distribute the state vector stack
                    rho=codistributed(rho,codist);
                    
                    % Get local part
                    rho_local=getLocalPart(rho);
                    
                    % Estimate the scaling coefficient
                    scaling=max([1 norm(rho_local,1)]);
                    
                    % Scale the vector and make it full
                    rho_local=full(rho_local/scaling);
                    
                    % Run the Krylov procedure
                    for n=1:nsteps
                        next_term=rho_local; k=1;
                        while nnz(abs(next_term)>eps)>0
                            next_term=-1i*(L*next_term)*time_step/(k*nsteps);
                            rho_local=rho_local+next_term; k=k+1;
                        end
                    end
                    
                    % Scale the vector back
                    rho_local=rho_local*scaling;
                
                end
                
                % Gather the array
                rho=[rho_local{:}];
                
            else
                
                % Estimate the scaling coefficient
                scaling=max([1 norm(rho,1)]);
                
                % Scale the vector and make it full
                rho=full(rho/scaling);
                
                % Move to GPU if needed
                if ismember('gpu',spin_system.sys.enable)
                    L=gpuArray(L); rho=gpuArray(rho);
                end
                
                % Run the Krylov procedure
                for n=1:nsteps
                    next_term=rho; k=1;
                    while nnz(abs(next_term)>eps)>0
                        next_term=-1i*(time_step/(k*nsteps))*(L*next_term);
                        rho=rho+next_term; k=k+1;
                    end
                end
                
                % Scale the vector back
                rho=gather(rho)*scaling;
                
            end
                
        end
        
    case 'zeeman-hilb'
        
        % Decide the propagator calculation method
        if size(L,1)>spin_system.tols.small_matrix
            
            % For large matrices and tensor trains use Spinach propagator module
            P=propagator(spin_system,L,time_step);
            
        else
            
            % For small matrices use Matlab's expm
            P=expm(-1i*L*time_step);
            
        end
        
        % Perform the propagation step
        if iscell(rho)
            parfor n=1:numel(rho)
                rho{n}=P*rho{n}*P';
            end
        else
            rho=P*rho*P';
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% Consistency enforcement
function grumble(L,time_step)
if (~isnumeric(L))||(~isnumeric(time_step))
    error('L and time_step must be numeric.');
end
if ~isscalar(time_step)
    error('timestep must be a scalar.');
end
end

% Evans boldly put 50 atm of ethylene in a cell with 25 atm of oxygen. The
% apparatus subsequently blew up, but luckily not before he had obtained
% the spectra shown in Figure 8.
%
% A.J. Mehrer and R.S. Mulliken, Chem. Rev. 69 (1969) 639-656.

