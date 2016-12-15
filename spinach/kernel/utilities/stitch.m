% Stitching function for bidirectionally propagated 3D experiments. Syntax:
%
% fid=stitch(spin_system,L,rho_stack,coil_stack,mtp_oper,mtp_time,t1,t2,t3)
%
% Parameters:
%
%             L - spin system Liouvillian
% 
%     rho_stack - state vector stack from the forward part of
%                 the simulation
%
%    coil_stack - coil vector stack from the backward part of
%                 the sumulation
%
%      mtp_oper - operator for the pulse to be executed in the
%                 middle of the t2 period
%
%      mtp_time - duration for the pulse to be executed in the
%                 middle of the t2 period
%
%     t1.nsteps - number of time steps in t1
%
%     t2.nsteps - number of time steps in t2
%
%     t3.nsteps - number of time steps in t3
%
% The function returns the three-dimensional free induction decay.
%
% i.kuprov@soton.ac.uk

function fid=stitch(spin_system,L,rho_stack,coil_stack,mtp_oper,mtp_time,t1,t2,t3)

% Preallocate the fid
fid=zeros(t3.nsteps,t2.nsteps,t1.nsteps);

% Run the reduction without splitting the space
spin_system.sys.disable={'pt','symmetry'};
P=reduce(spin_system,L+mtp_oper,[rho_stack coil_stack]);
rho_stack=P{1}'*rho_stack; coil_stack=P{1}'*coil_stack;
L=P{1}'*L*P{1}; mtp_oper=P{1}'*mtp_oper*P{1};

% Compute propagators
P_forw=propagator(spin_system,L,t2.timestep/2); P_back=P_forw';
P_mtp=propagator(spin_system,mtp_oper,mtp_time);

% Stitch the trajectories
for k=1:t2.nsteps
    report(spin_system,['stitching forward and backward trajectories, step '...
           num2str(k) '/' num2str(t2.nsteps) '...']);
    fid(:,k,:)=coil_stack'*P_mtp*rho_stack;
    rho_stack=P_forw*rho_stack; 
    coil_stack=P_back*coil_stack;
end

end

% A wolf hates both men and dogs, but dogs he hates more.
%
% Sergey Dovlatov

