% Krotov type algorithms for phase-sensitive state-to-state transfer prob-
% lems.
%
% <http://spindynamics.org/wiki/index.php?title=Krotov.m>

function [diag_data,waveform,traj_fwd,traj_bwd,fidelity,grad]=...
    krotov(spin_system,ctrl_system,drift,waveform,traj_fwd,traj_bwd,hess)

% Scale the waveform
waveform=ctrl_system.power_level*waveform;

% scale the reference waveform (if supplied) else use last waveform
if ~isfield(ctrl_system,'ref_waveform')
    ctrl_system.ref_waveform=waveform;
end

% Silence spinach output for propagator diagnostics
outtemp=spin_system.sys.output;
spin_system.sys.output='hush';

if exist('hess','var')
    hess_index=reshape(1:numel(waveform),[numel(ctrl_system.control_ops) ctrl_system.nsteps]);
end

if isempty(traj_fwd)
    
    % initialise trajectories
    traj_fwd=zeros(size(ctrl_system.rho,1),ctrl_system.nsteps+1);
    traj_bwd=zeros(size(ctrl_system.rho,1),ctrl_system.nsteps+1);
    
    % Run the initial forward propagation
    traj_fwd(:,1)=ctrl_system.rho;
    for n=1:ctrl_system.nsteps
        L=drift;
        for k=1:numel(ctrl_system.control_ops)
            L=L+waveform(k,n)*ctrl_system.control_ops{k};
        end
        traj_fwd(:,n+1)=step(spin_system,L,traj_fwd(:,n),ctrl_system.time_step);
    end
    
    % Run the initial backward propagation
    traj_bwd(:,ctrl_system.nsteps+1)=ctrl_system.target;
    for n=ctrl_system.nsteps:-1:1
        L=drift';
        for k=1:numel(ctrl_system.control_ops)
            L=L+waveform(k,n)*ctrl_system.control_ops{k};
        end
        traj_bwd(:,n)=step(spin_system,L,traj_bwd(:,n+1),-ctrl_system.time_step);
    end
    
    
else
    
    % Run forward Krotov sweep
    traj_fwd(:,1)=ctrl_system.rho;
    for n=1:ctrl_system.nsteps
        
        % Run microiterations
        difference=1;
        while difference > ctrl_system.krotov_microiteration
            
            % Start the Liouvillian and store previous control values
            L=drift; prev_values=waveform(:,n);
            
            for k=1:numel(ctrl_system.control_ops)
                if exist('hess','var')
                    waveform(k,n)=ctrl_system.ref_waveform(k,n)+...
                        (1/ctrl_system.lambda)*(1/hess(hess_index(k,n),hess_index(k,n)))*imag(traj_bwd(:,n+1)'*ctrl_system.control_ops{k}*traj_fwd(:,n));
                    for p=(n-1):-1:1
                        waveform(k,n)=waveform(k,n)-...
                            (1/hess(hess_index(k,n),hess_index(k,n)))*hess(hess_index(k,n),hess_index(k,p))*(waveform(k,p)-ctrl_system.ref_waveform(k,p));
                    end
                else
                    waveform(k,n)=ctrl_system.ref_waveform(k,n)+...
                        (1/ctrl_system.lambda)*imag(traj_bwd(:,n+1)'*ctrl_system.control_ops{k}*traj_fwd(:,n));
                end
                L=L+waveform(k,n)*ctrl_system.control_ops{k};
            end
            
            % Compute the difference
            difference=norm(waveform(:,n)-prev_values);
            
            % Refresh the forward point
            traj_fwd(:,n+1)=step(spin_system,L,traj_fwd(:,n),ctrl_system.time_step);
        end
    end
    
    % Run backward Krotov sweep
    traj_bwd(:,ctrl_system.nsteps+1)=ctrl_system.target;
    
    for n=ctrl_system.nsteps:-1:1
        L=drift';
        for k=1:numel(ctrl_system.control_ops)
            if ismember(ctrl_system.method,{'zhu-rabitz'})
                waveform(k,n)=waveform(k,n)+(1/ctrl_system.lambda)*imag(traj_bwd(:,n+1)'*ctrl_system.control_ops{k}*traj_fwd(:,n+1));
            end
            L=L+waveform(k,n)*ctrl_system.control_ops{k};
        end
        traj_bwd(:,n)=step(spin_system,L,traj_bwd(:,n+1),-ctrl_system.time_step);
    end
    
end

% compute the gradient if needed
if exist('hess','var')
    % Preallocate results
    grad=zeros(size(waveform));
    for n=1:ctrl_system.nsteps
        for k=1:numel(ctrl_system.control_ops)
            grad(k,n)=-2*imag(traj_bwd(:,n+1)'*ctrl_system.control_ops{k}*traj_fwd(:,n))-...
                2*ctrl_system.power_level*ctrl_system.lambda*(waveform(k,n)-ctrl_system.ref_waveform(k,n));
        end
    end
    grad=grad./ctrl_system.power_level;
end

% Normalize the total objective, total gradient, and total Hessian
fidelity=real(ctrl_system.target'*traj_fwd(:,end));
waveform=real(waveform./ctrl_system.power_level);

% re-enable spinach console output
spin_system.sys.output=outtemp;

% Compile diagnostics data
if nargout>1
    diag_data.rho=ctrl_system.rho;
    diag_data.target=ctrl_system.target;
    diag_data.trajectory=traj_fwd;
    diag_data.trajectory_bwd=traj_bwd;
    diag_data.drift=drift;
end

end

% It wasn't a dark and stormy night. It should have been, but there's
% the weather for you. For every mad scientist who's had a convenient
% thunderstorm just on the night his Great Work is complete and lying
% on the slab, there have been dozens who've sat around aimlessly un-
% der the peaceful stars while Igor clocks up the overtime.
%
% Terry Pratchett