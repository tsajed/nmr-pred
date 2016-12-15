% Shaped pulse in amplitude-frequency coordinates using Fokker-Planck
% formalism. Syntax:
%
%   rho=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,...
%                   rf_amp_list,rf_dur_list,rf_phi,max_rank,method)
%
% Parameters:
%
%        L0          - drift Liouvillian that continues
%                      running in the background
%
%        Lx          - X projection of the RF operator
%
%        Ly          - Y projection of the RF operator
%
%        rho         - initial condition
%
%        rf_frq_list - a vector of RF frequencies at each
%                      time slice, Hz
%
%        rf_amp_list - a vector of RF amplitudes at each
%                      time slice, Hz
%
%        rf_dur_list - a vector of time slice durations,
%                      in seconds
%
%        rf_phi      - RF phase at time zero
%
%        max_rank    - maximum rank of the Fokker-Planck
%                      theory
%
%        method      - propagation method, 'expv' for Krylov
%                      propagation, 'expm' for exponential
%                      propagation, 'evolution' for Spinach
%                      evolution function
%
% i.kuprov@soton.ac.uk

function rho=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,rf_dur_list,rf_phi,max_rank,method)

% Set the defaults
if ~exist('method','var'), method='expv'; end

% Get problem dimensions
spc_dim=2*max_rank+1; spn_dim=size(L0,1); stk_dim=size(rho,2);
report(spin_system,['lab space problem dimension     ' num2str(spc_dim)]);
report(spin_system,['spin space problem dimension    ' num2str(spn_dim)]);
report(spin_system,['state vector stack size         ' num2str(stk_dim)]);
report(spin_system,['Fokker-Planck problem dimension ' num2str(spc_dim*spn_dim)]);

% Compute RF angles and Fourier derivative operator
[rotor_angles,d_dphi]=fourdif(spc_dim,1);

% Add the phase
rotor_angles=rotor_angles+rf_phi;

% Build the background Liouvillian
F0=kron(speye(spc_dim),L0);

% Build the RF operator
F1=cell(spc_dim);
for n=1:spc_dim
    for k=1:spc_dim
        F1{n,k}=spalloc(spn_dim,spn_dim,0);
    end
    F1{n,n}=cos(rotor_angles(n))*Lx+sin(rotor_angles(n))*Ly;
end
F1=clean_up(spin_system,cell2mat(F1),spin_system.tols.liouv_zero);

% Build the phase increment operator
M=kron(d_dphi,speye(size(L0)));

% Project the state
P=spalloc(spc_dim,1,1); P(1)=1; rho=kron(P,rho);

% Run the pulse
switch method
    
    case 'expv'
        
        % Use Krylov propagation
        for n=1:numel(rf_frq_list)
            report(spin_system,['Krylov propagation step ' num2str(n) ' of ' num2str(numel(rf_frq_list)) '...']);
            rho=step(spin_system,F0+2*pi*rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,rho,rf_dur_list(n));
        end
        
    case 'expm'
        
        % Use exponential propagation
        for n=1:numel(rf_frq_list)
            rho=propagator(spin_system,F0+2*pi*rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,rf_dur_list(n))*rho;
        end
        
    case 'evolution'
        
        % Use the evolution function
        for n=1:numel(rf_frq_list)
            rho=evolution(spin_system,F0+2*pi*rf_amp_list(n)*F1+2i*pi*rf_frq_list(n)*M,[],rho,rf_dur_list(n),1,'final');
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown propagation method.');
        
end

% Fold back the state
rho=squeeze(sum(reshape(full(rho),[spn_dim spc_dim stk_dim]),2));

end

% A man gazing at the stars is at the mercy 
% of the puddles in the road.
%
% Alexander Smith

