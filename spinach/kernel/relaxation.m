% The relaxation superoperator. All options are set during the call to
% create.m function. Syntax:
%
%                 R=relaxation(spin_system,euler_angles)
%
% where euler_angles should be used for solid state systems in which
% the thermal equilibrium state is orientation-dependent.
%
% Further information is available in the following papers:
%
%         http://dx.doi.org/10.1016/j.jmr.2010.12.004 (Redfield)
%         http://dx.doi.org/10.1002/mrc.1242          (Lindblad)
%         http://dx.doi.org/10.1007/s00723-012-0367-0 (Nottingham DNP)
%         http://dx.doi.org/10.1039/c2cp42897k        (Weizmann DNP)
%
% Note: Nottingham relaxation theory is only applicable to cross effect
%       DNP systems and Weizmann relaxation theory to solid effect and
%       cross effect DNP systems. 
%
% Note: The function has been written for minimal memory footprint and
% performs aggressive memory recycling. Relaxation superoperator dimen-
% sions in excess of 1,000,000 are in practice feasible.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function R=relaxation(spin_system,euler_angles)

% Set defaults
spin_system=defaults(spin_system);

% Check consistency
grumble(spin_system);

% Get the matrix going
R=mprealloc(spin_system,0);

% Do nothing if none specified
if isempty(spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'relaxation superoperator set to zero.');
    
    % Exit the function
    return;
    
end

% Add damping terms
if ismember('damp',spin_system.rlx.theories)
        
    % Set up damping
    switch spin_system.bas.formalism
        
        case {'zeeman-hilb'}
            
            % Update the user
            report(spin_system,['isotropic damping at ' num2str(spin_system.rlx.damp_rate) ' Hz for all states.']);
            
            % Damp everything at the rate specified
            R=R-(spin_system.rlx.damp_rate/2)*unit_oper(spin_system);
            
        case {'zeeman-liouv','sphten-liouv'}
            
            % Update the user
            report(spin_system,['isotropic damping at ' num2str(spin_system.rlx.damp_rate) ' Hz for all states except the unit state.']);
            
            % Damp everything at the rate specified
            R=R-spin_system.rlx.damp_rate*unit_oper(spin_system);
            
            % Make sure the unit state is not damped
            U=unit_state(spin_system);
            R=R-(U'*R*U)*(U*U');
            
    end
    
    % Inform the user
    report(spin_system,'isotropic damping terms have been incorporated.');
    
end

% Add H-strain terms
if ismember('hstrain',spin_system.rlx.theories)
    
    % Compute the damp rate for the current orientation
    n=euler2dcm(euler_angles(1),euler_angles(2),euler_angles(3))*[0; 0; 1];
    hstrain_rates=sqrt(sum((spin_system.rlx.hstrain_rates.^2)*(n.^2)));
    
    % Set up damping
    switch spin_system.bas.formalism
        
        case {'zeeman-hilb'}
            
            % Update the user
            report(spin_system,['anisotropic damping at ' num2str(spin_system.rlx.damp_rate) ' Hz for all states.']);
            
            % Damp everything at the rate specified
            R=R-(hstrain_rates/2)*unit_oper(spin_system);
            
        case {'zeeman-liouv','sphten-liouv'}
            
            % Update the user
            report(spin_system,['anisotropic damping at ' num2str(spin_system.rlx.damp_rate) ' Hz for all states except the unit state.']);
            
            % Damp everything at the rate specified
            R=R-hstrain_rates*unit_oper(spin_system);
            
            % Make sure the unit state is not damped
            U=unit_state(spin_system);
            R=R-(U'*R*U)*(U*U');
            
    end
    
    % Inform the user
    report(spin_system,'anisotropic damping terms have been incorporated.');
    
end
    
% Add T1/T2 terms
if ismember('t1_t2',spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'extended T1,T2 relaxation theory with user-supplied R1 and R2 rates.');
    
    % Prealocate the relaxation rate array
    matrix_dim=size(spin_system.bas.basis,1);
    relaxation_rates=zeros(matrix_dim,1);
    
    % Compute the ranks and projections
    [L,M]=lin2lm(spin_system.bas.basis);
    
    % Extract rates
    r1_rates=spin_system.rlx.r1_rates;
    r2_rates=spin_system.rlx.r2_rates;
    
    % Inspect every state and assign the relaxation rate
    parfor n=1:matrix_dim
        
        % Copy rate vectors to nodes
        local_r1_rates=r1_rates;
        local_r2_rates=r2_rates;
        
        % Spins in unit state do not contribute
        mask=(L(n,:)~=0);
        
        % Spins in longitudinal states contribute their R1
        r1_spins=(~logical(M(n,:)))&mask;
        r1_sum=sum(local_r1_rates(r1_spins));
        
        % Spins in transverse states contribute their R2
        r2_spins=(logical(M(n,:)))&mask;
        r2_sum=sum(local_r2_rates(r2_spins));
        
        % Total relaxation rate for the state
        relaxation_rates(n)=r1_sum+r2_sum;
        
    end
    
    % Deallocate the rank and projection arrays
    clear('L','M');
    
    % Form the relaxation superoperator
    R=R-spdiags(relaxation_rates,0,matrix_dim,matrix_dim);
    
    % Inform the user
    report(spin_system,'extended T1,T2 terms have been incorporated.');
    
end

% Add Redfield terms
if ismember('redfield',spin_system.rlx.theories)&&(norm(spin_system.rlx.tau_c)>0)
    
    % Update the user
    report(spin_system,'Redfield theory with user-supplied anisotropies and correlation times.');
    
    % Get the rotational basis, including the non-secular terms
    report(spin_system,'computing the lab frame Hamiltonian superoperator...');
    [L0,Q]=hamiltonian(assume(spin_system,'labframe'));
    
    % Kill the terms in the static Liouvillian that are irrelevant on the time scale of the integration
    timescale=max(spin_system.rlx.tau_c)*log(1/spin_system.tols.rlx_integration);
    L0=clean_up(spin_system,L0,spin_system.tols.rlx_integration/timescale);
    
    % Run the calculation in parallel
    report(spin_system,'codistributing Hamiltonian components...');
    spmd
        
        % Distribute the work
        kmpq=codistributed(1:625,codistributor1d(2));
        
        % Preallocate local answer
        R_loc=mprealloc(spin_system,10);
        
        % Loop over the local multi-index array
        for kmpq_local=getLocalPart(kmpq)
            
            % Extract indices
            [k,m,p,q]=ind2sub([5 5 5 5],kmpq_local);
            
            % Get decay weights and rates for correlation function
            [weights,rates]=corrfun(spin_system,k,m,p,q);
            
            % Loop over correlation function exponentials
            for j=1:numel(weights)
                
                % Compute the term in the rotational expansion sum
                if significant(Q{k,m},eps)&&significant(Q{p,q},eps)&&significant(weights(j),eps)
                    
                    % Inform the user
                    report(spin_system,['integrating spherical component k=' num2str(3-k,'%+d') ', m=' num2str(3-m,'%+d') ...
                                                                      ', p=' num2str(3-p,'%+d') ', q=' num2str(3-q,'%+d') '...']);
                
                    % Set the upper integration limit according to the accuracy goal
                    upper_limit=-1.5*(1/rates(j))*log(1/spin_system.tols.rlx_integration);
                                                                      
                    % Compute the integral using the auxiliary matrix method
                    R_loc=R_loc-weights(j)*Q{k,m}*expmint(spin_system,L0,Q{p,q}',L0-1i*rates(j)*speye(size(L0)),upper_limit);
                        
                end
                
            end
            
        end
        
    end
    
    % Deallocate variables
    clear('Q','L0');
    
    % Gather data from the nodes
    report(spin_system,'gathering data from the nodes...');
    for n=1:numel(R_loc), R=R+R_loc{n}; R_loc{n}=[]; end
    
    % Update the user
    report(spin_system,'Redfield theory terms have been incorporated.');
    
end
        
% Add Lindblad terms
if ismember('lindblad',spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'Lindblad theory with user-supplied R1 and R2 rates.');
    
    % Loop over spins
    for n=1:spin_system.comp.nspins
        
        % Get the basic superoperators
        Lp_left=operator(spin_system,{'L+'},{n},'left');
        Lp_right=operator(spin_system,{'L+'},{n},'right');
        Lm_left=operator(spin_system,{'L-'},{n},'left');
        Lm_right=operator(spin_system,{'L-'},{n},'right');
        Lz_left=operator(spin_system,{'Lz'},{n},'left');
        Lz_right=operator(spin_system,{'Lz'},{n},'right');
        
        % Get the rates
        Rpm=spin_system.rlx.lind_r1_rates(n)/4;
        Rmp=spin_system.rlx.lind_r1_rates(n)/4;
        Rzz=spin_system.rlx.lind_r2_rates(n)-...
            spin_system.rlx.lind_r1_rates(n)/2;
        
        % Build the superoperator
        R=R+Rpm*(2*Lp_left*Lm_right-Lm_left*Lp_left-Lp_right*Lm_right)+...
            Rmp*(2*Lm_left*Lp_right-Lp_left*Lm_left-Lm_right*Lp_right)+...
            Rzz*(2*Lz_left*Lz_right-Lz_left*Lz_left-Lz_right*Lz_right);
        
    end
    
    % Inform the user
    report(spin_system,'Lindblad theory terms have been incorporated.');
    
end
    
% Add Weizmann DNP theory terms
if ismember('weizmann',spin_system.rlx.theories)
    
    % Update the user
    report(spin_system,'Weizmann DNP theory with user-supplied R1e, R2e, R1n, R2n, R1d and R2d rates.');
    
    % Get electron superoperators
    Ep_L=operator(spin_system,'L+','electrons','left');
    Ep_R=operator(spin_system,'L+','electrons','right');
    Em_L=operator(spin_system,'L-','electrons','left');
    Em_R=operator(spin_system,'L-','electrons','right');
    Ez_L=operator(spin_system,'Lz','electrons','left');
    Ez_R=operator(spin_system,'Lz','electrons','right');
    
    % Get electron rates
    Rpm=spin_system.rlx.weiz_r1e/4; 
    Rmp=spin_system.rlx.weiz_r1e/4;
    Rzz=spin_system.rlx.weiz_r2e-...
        spin_system.rlx.weiz_r1e/2;
    
    % Add Lindblad type terms
    R=R+Rpm*(2*Ep_L*Em_R-Em_L*Ep_L-Ep_R*Em_R)+...
        Rmp*(2*Em_L*Ep_R-Ep_L*Em_L-Em_R*Ep_R)+...
        Rzz*(2*Ez_L*Ez_R-Ez_L*Ez_L-Ez_R*Ez_R);
    
    % Get nuclear superoperators
    Np_L=operator(spin_system,'L+','nuclei','left');
    Np_R=operator(spin_system,'L+','nuclei','right');
    Nm_L=operator(spin_system,'L-','nuclei','left');
    Nm_R=operator(spin_system,'L-','nuclei','right');
    Nz_L=operator(spin_system,'Lz','nuclei','left');
    Nz_R=operator(spin_system,'Lz','nuclei','right');
    
    % Get nuclear rates
    Rpm=spin_system.rlx.weiz_r1n/4;
    Rmp=spin_system.rlx.weiz_r1n/4;
    Rzz=spin_system.rlx.weiz_r2n-...
        spin_system.rlx.weiz_r1n/2;
    
    % Add Lindblad type terms
    R=R+Rpm*(2*Np_L*Nm_R-Nm_L*Np_L-Np_R*Nm_R)+...
        Rmp*(2*Nm_L*Np_R-Np_L*Nm_L-Nm_R*Np_R)+...
        Rzz*(2*Nz_L*Nz_R-Nz_L*Nz_L-Nz_R*Nz_R);
    
    % Add dipolar cross-relaxation terms
    for n=1:spin_system.comp.nspins
        for k=1:spin_system.comp.nspins
            
            % Process longitudinal terms
            if spin_system.rlx.weiz_r1d(n,k)~=0
                
                % Generate flip-flop superoperator
                FF_L=operator(spin_system,{'L+','L-'},{n,k},'left')+...
                     operator(spin_system,{'L-','L+'},{n,k},'left');
                FF_R=operator(spin_system,{'L+','L-'},{n,k},'right')+...
                     operator(spin_system,{'L-','L+'},{n,k},'right');
                
                % Add Lindblad type term
                R=R+spin_system.rlx.weiz_r1d(n,k)*(2*FF_L*FF_R-FF_L*FF_L-FF_R*FF_R)/2;
                
            end
            
            % Process transverse terms
            if spin_system.rlx.weiz_r2d(n,k)~=0
                
                % Generate ZZ superoperator
                ZZ_L=operator(spin_system,{'Lz','Lz'},{n,k},'left');
                ZZ_R=operator(spin_system,{'Lz','Lz'},{n,k},'right');
                
                % Add Lindblad type term
                R=R+spin_system.rlx.weiz_r2d(n,k)*(2*ZZ_L*ZZ_R-ZZ_L*ZZ_L-ZZ_R*ZZ_R)/2;
                
            end
            
        end
    end
    
    % Inform the user
    report(spin_system,'Weizmann DNP relaxation theory terms have been incorporated.');
    
end
        
% Add Nottingham DNP theory terms
if ismember('nottingham',spin_system.rlx.theories)
        
    % Update the user
    report(spin_system,'Nottingham DNP theory with user-supplied R1e, R1n, R2e and R2n rates.');
    
    % Set shorthands
    r1e=spin_system.rlx.nott_r1e; r2e=spin_system.rlx.nott_r2e;
    r1n=spin_system.rlx.nott_r1n; r2n=spin_system.rlx.nott_r2n;
    
    % Find electrons
    electrons=find(strcmp('E',spin_system.comp.isotopes));
    
    % Get basic electron operators
    S1p=operator(spin_system,{'L+'},{electrons(1)});
    S1z=operator(spin_system,{'Lz'},{electrons(1)});
    S2p=operator(spin_system,{'L+'},{electrons(2)});
    S2z=operator(spin_system,{'Lz'},{electrons(1)});
    S1zS2z=operator(spin_system,{'Lz','Lz'},{electrons(1),electrons(2)});
    S1zS2p=operator(spin_system,{'Lz','L+'},{electrons(1),electrons(2)});
    S1pS2z=operator(spin_system,{'L+','Lz'},{electrons(1),electrons(2)});
    
    % Form component operators
    O11=0.5*(+S1z+S2z)+S1zS2z; O12=0.5*S2p+S1zS2p;
    O22=0.5*(+S1z-S2z)-S1zS2z; O13=0.5*S1p+S1pS2z;
    O33=0.5*(-S1z+S2z)-S1zS2z; O24=0.5*S1p-S1pS2z;
    O44=0.5*(-S1z-S2z)+S1zS2z; O34=0.5*S2p-S1zS2p;
        
    % Build the electron part of the relaxation superoperator
    R=-0.25*r1e*(O12*O12' + O12'*O12 - O11*O11 - O22*O22 + O13*O13' + O13'*O13 - O11*O11 - O33*O33 +...
                 O24*O24' + O24'*O24 - O22*O22 - O44*O44 + O34*O34' + O34'*O34 - O33*O33 - O44*O44)+...
       0.5*(r2e*(O11*O22  + O22*O11  + O11*O33 + O33*O11 + O44*O22  + O22*O44  + O33*O44 + O44*O33)+...
            r2n*(O11*O44  + O44*O11  + O22*O33 + O33*O22));
        
    % Find nuclei
    nuclei=find(~cellfun(@(x)strncmp(x,'E',1),spin_system.comp.isotopes));
        
    % Loop over nuclei
    for n=nuclei
            
        % Get basic nuclear operators
        Iz=operator(spin_system,{'Lz'},{n});
        Ip=operator(spin_system,{'L+'},{n});
        
        % Build the nuclear part of the relaxation superoperator
        R=R-0.25*r1n*(Ip*Ip'+Ip'*Ip)-r2n*Iz*Iz;
        
    end
    
    % Inform the user
    report(spin_system,'Nottingham DNP relaxation theory terms have been incorporated.');
        
end

% Add scalar relaxation of the first kind
if ismember('SRFK',spin_system.rlx.theories)&&...
           (spin_system.rlx.srfk_tau_c>0)

    % Inform the user
    report(spin_system,'scalar relaxation of the first kind...');
    
    % Set local assumptions
    spin_system_local=assume(spin_system,spin_system.rlx.srfk_assume);
    
    % Get the background Hamiltonian
    H0=hamiltonian(spin_system_local);
    
    % Set local assumptions
    spin_system_local=assume(spin_system,spin_system.rlx.srfk_assume,'couplings');
    
    % Replace couplings with their modulation depths
    [rows,cols]=find(~cellfun(@isempty,spin_system.rlx.srfk_mdepth));
    spin_system_local.inter.coupling.matrix=cell(size(spin_system.inter.coupling.matrix));
    for n=1:numel(rows)
        spin_system_local.inter.coupling.matrix{rows(n),cols(n)}=...
             2*pi*eye(3)*spin_system.rlx.srfk_mdepth{rows(n),cols(n)};
    end
    
    % Get the modulated coupling Hamiltonian
    H1=hamiltonian(spin_system_local);
    
    % Set the upper integration limit according to the accuracy goal
    upper_limit=2*spin_system.rlx.srfk_tau_c*log(1/spin_system.tols.rlx_integration);
    
    % Kill the terms in H0 that are irrelevant on the time scale of the integration
    H0=clean_up(spin_system,H0,1e-2/upper_limit);
    
    % Take the integral using the auxiliary matrix exponential technique
    report(spin_system,'integrating the SRFK component...');
    decay_rate=-1/spin_system.rlx.srfk_tau_c;
    R=R-H1*expmint(spin_system,H0,H1',H0-1i*decay_rate*speye(size(H0)),upper_limit);
    
    % Deallocate variables
    clear('H0','H1','spin_system_local');
    
    % Inform the user
    report(spin_system,'SRFK terms have been incorporated.');
    
end

% Add scalar relaxation of the second kind
if ismember('SRSK',spin_system.rlx.theories)
    
    % Inform the user
    report(spin_system,'scalar relaxation of the second kind...');
    
    % Preallocate the answers
    R1=zeros(spin_system.comp.nspins); R2=zeros(spin_system.comp.nspins);
    
    % Loop over pairs of spins
    for n=1:spin_system.comp.nspins
        for k=1:spin_system.comp.nspins
            
            % Determine the J-coupling
            J=0;
            if ~isempty(spin_system.inter.coupling.matrix{n,k})
                J=J+trace(spin_system.inter.coupling.matrix{n,k})/3;
            end
            if ~isempty(spin_system.inter.coupling.matrix{k,n})
                J=J+trace(spin_system.inter.coupling.matrix{k,n})/3;
            end

            % Only process heteronuclear pairs
            if (~strcmp(spin_system.comp.isotopes{n},spin_system.comp.isotopes{k}))&&(n>k)&&(abs(J)>0)
                
                % Find relaxation rates for the spins in question
                Lz_n=state(spin_system,{'Lz'},{n}); Lz_n=Lz_n/norm(Lz_n);
                Lp_n=state(spin_system,{'L+'},{n}); Lp_n=Lp_n/norm(Lp_n);
                R1n=-Lz_n'*R*Lz_n; R2n=-Lp_n'*R*Lp_n;
                Lz_k=state(spin_system,{'Lz'},{k}); Lz_k=Lz_k/norm(Lz_k);
                Lp_k=state(spin_system,{'L+'},{k}); Lp_k=Lp_k/norm(Lp_k);
                R1k=-Lz_k'*R*Lz_k; R2k=-Lp_k'*R*Lp_k;
                
                % Find multiplicities
                In=spin_system.comp.mults(n); Ik=spin_system.comp.mults(k);
                
                % Find frequency difference
                delta_omega=spin_system.inter.basefrqs(n)-spin_system.inter.basefrqs(k);
        
                % Compute Abragam's expressions
                R1(n)=R1(n)+(8/3)*(J^2)*Ik*(Ik+1)*(R2k/(R2k^2+delta_omega^2));
                R1(k)=R1(k)+(8/3)*(J^2)*In*(In+1)*(R2n/(R2n^2+delta_omega^2));
                R2(n)=R2(n)+(4/3)*(J^2)*Ik*(Ik+1)*(1/R1k+R2k/(R2k^2+delta_omega^2));
                R2(k)=R2(k)+(4/3)*(J^2)*Ik*(In+1)*(1/R1n+R2n/(R2n^2+delta_omega^2));
                
            end
            
        end 
    end
    
    % Clone the spin system structure
    spin_system_local=spin_system;
    spin_system_local.rlx.theories={'t1_t2'};
    spin_system_local.rlx.r1_rates=R1;
    spin_system_local.rlx.r2_rates=R2;
    
    % Issue a recursive call for SRSK
    report(spin_system,'recursive call for SRSK relaxation terms...'); 
    R=R+relaxation(spin_system_local);
    
    % Inform the user
    report(spin_system,'SRSK relaxation theory terms have been incorporated.');
    
end

% Print matrix density statistics
report(spin_system,['full relaxation superoperator density ' num2str(100*nnz(R)/numel(R))...
                    '%, nnz ' num2str(nnz(R)) ', sparsity ' num2str(issparse(R))]);

% Decide the fate of dynamic frequency shifts
switch spin_system.rlx.dfs
    
    case 'keep'
        
        % Do nothing and inform the user
        report(spin_system,'dynamic frequency shifts have been kept.');
        
    case 'ignore'
        
        % Kill the dynamic frequency shifts
        R=real(R);
        
        % Inform the user
        report(spin_system,'dynamic frequency shifts have been ignored.');
        
    otherwise
        
        % Complain and bomb out
        error('invalid value of the inter.rlx_dfs parameter.');
        
end

% Decide the fate of the off-diagonal components
switch spin_system.rlx.keep
    
    case 'diagonal'
        
        % Pull out the diagonal
        R=diag(diag(R));
        
        % Inform the user
        report(spin_system,'all cross-relaxation terms have been ignored.');
        
    case 'kite'
        
        % Refuse to process inappropriate cases
        if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
            error('kite option is only available for sphten-liouv formalism.');
        end
        
        % Compile the index of all longitudinal spin orders
        [~,M]=lin2lm(spin_system.bas.basis);
        long_states=find(sum(abs(M),2)==0);
        
        % Index the relaxation superoperator
        [rows,cols,vals]=find(R);
        
        % Keep self-relaxation and longitudinal cross-relaxation terms
        vals=vals.*((ismember(rows,long_states)&ismember(cols,long_states))|(rows==cols));
        
        % Recompose the relaxation superoperator
        R=sparse(rows,cols,vals,length(R),length(R)); clear('rows','cols','vals');
        
        % Inform the user
        report(spin_system,'transverse cross-relaxation terms have been ignored.');
        
    case 'secular'
        
        % Refuse to process inappropriate cases
        if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
            error('secular option is only available for sphten-liouv formalism.');
        end
        
        % Compute state carrier frequencies
        [~,M]=lin2lm(spin_system.bas.basis);
        matrix_dim=size(spin_system.bas.basis,1);
        frequencies=sum(repmat(spin_system.inter.basefrqs,matrix_dim,1).*M,2);
        
        % Index the relaxation superoperator
        [rows,cols,vals]=find(R);
        
        % Set the initial keep mask
        keep_mask=false(size(vals));
        
        % Loop over unique frequencies
        for omega=unique(frequencies)'
            
            % Find the states having the current frequency
            current_frq_group=find(frequencies==omega);
            
            % Update the keep mask
            keep_mask=keep_mask|(ismember(rows,current_frq_group)&...
                                 ismember(cols,current_frq_group));
            
        end
        
        % Recompose the relaxation superoperator
        R=sparse(rows,cols,vals.*keep_mask,length(R),length(R)); clear('rows','cols','vals');
        
        % Inform the user
        report(spin_system,'non-secular cross-relaxation terms have been ignored.');
        
    case 'labframe'
        
        % Inform the user
        report(spin_system,'returning complete relaxation superoperator (lab frame simulations only).');
        
end
    
% Print matrix density statistics
report(spin_system,['final relaxation superoperator density ' num2str(100*nnz(R)/numel(R))...
                    '%, nnz(R)=' num2str(nnz(R)) ', sparsity ' num2str(issparse(R))]);
    
% Choose the thermalization model
switch spin_system.rlx.equilibrium
    
    case 'zero'
        
        % Do nothing and inform the user
        report(spin_system,'WARNING - the spin system will relax to the all-zero state.');
        
    case 'levante'
        
        % Inform the user
        report(spin_system,'the system will relax to thermal equilibrium (Levante-Ernst formalism).');
        
        % Use Levante-Ernst thermalization
        if exist('euler_angles','var')
            [H,Q]=hamiltonian(assume(spin_system,'labframe'),'left');
            R(:,1)=-R*equilibrium(spin_system,H,Q,euler_angles);
            report(spin_system,'equilibrium state was computed using anisotropic Hamiltonian at');
            report(spin_system,['alpha=' num2str(euler_angles(1)) ', beta=' num2str(euler_angles(2)) ', gamma=' num2str(euler_angles(3))]);
        else
            H=hamiltonian(assume(spin_system,'labframe'),'left');
            R(:,1)=-R*equilibrium(spin_system,H);
            report(spin_system,'equilibrium state was computed using isotropic Hamiltonian.');
        end
        
    case 'dibari'
        
        % Inform the user
        report(spin_system,'the system will relax to thermal equilibrium (DiBari-Levitt formalism).');
        
        % Get the temperature factor
        beta=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature);
        
        % Use DiBari-Levitt thermalization
        if exist('euler_angles','var')
            [H,Q]=hamiltonian(assume(spin_system,'labframe'),'right');
            R=R*propagator(spin_system,H+orientation(Q,euler_angles),1i*beta);
            report(spin_system,'equilibrium state was computed using anisotropic Hamiltonian at');
            report(spin_system,['alpha=' num2str(euler_angles(1)) ', beta=' num2str(euler_angles(2)) ', gamma=' num2str(euler_angles(3))]);
        else
            H=hamiltonian(assume(spin_system,'labframe'),'right');
            R=R*propagator(spin_system,H,1i*beta);
            report(spin_system,'equilibrium state was computed using isotropic Hamiltonian.');
        end
        
    otherwise
        
        % Complain and bomb out
        error('unknown equilibrium specification.');
        
end

% Perform final clean-up
R=clean_up(spin_system,R,spin_system.tols.liouv_zero);

end

% Default option values
function spin_system=defaults(spin_system)
if ~isfield(spin_system.rlx,'equilibrium')
    report(spin_system,'relaxation destination not specified, assuming zero.');
    spin_system.rlx.equilibrium='zero';
end
if ~isfield(spin_system.rlx,'dfs')
    report(spin_system,'dynamic frequency shift policy not specified, DFS will be ignored.');
    spin_system.rlx.dfs='ignore';
end
end

% Consistency enforcement
function grumble(spin_system)
if ~isfield(spin_system,'rlx')
    error('relaxation data (.rlx) is missing from the spin_system structure.');
end
if ~isfield(spin_system.rlx,'theories')
    error('relaxation data (.rlx.theories) is missing from the spin_system structure.');
end
if ~iscell(spin_system.rlx.theories)||any(~cellfun(@ischar,spin_system.rlx.theories))
    error('spin_system.rlx.theories must be a cell array of character strings.');
end
for n=1:numel(spin_system.rlx.theories)
    if ~ismember(spin_system.rlx.theories{n},{'none','damp','hstrain','t1_t2','SRSK','SRFK'...
                                              'redfield','lindblad','nottingham','weizmann'})
        error('unrecognised relaxation theory specification.');
    end
end
if ( ismember('t1_t2',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv'}))
    error('extended T1,T2 relaxation theory is only available for sphten-liouv formalism.');
end
if ( ismember('redfield',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Redfield relaxation theory is only available in Liouville space.');
end
if ( ismember('lindblad',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Lindblad relaxation theory is only available in Liouville space.');
end
if ( ismember('nottingham',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Nottingham relaxation theory is only available in Liouville space.');
end
if ( ismember('weizmann',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('Weizmann relaxation theory is only available in Liouville space.');
end
if ( ismember('SRFK',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('SRFK relaxation theory is only available in Liouville space.');
end
if ( ismember('SRSK',spin_system.rlx.theories))&&...
   (~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'}))
    error('SRSK relaxation theory is only available in Liouville space.');
end
if ( strcmp(spin_system.rlx.equilibrium,'levante'))&&...
   (~strcmp(spin_system.bas.formalism,'sphten-liouv'))
    error('Levante-Ernst thermalization is only available for sphten-liouv formalism.');
end
if (strcmp(spin_system.rlx.equilibrium,'dibari'))&&...
   (spin_system.rlx.temperature==0)
    error('DiBari-Levitt thermalization cannot be used with the high-temperature approximation.')
end
end

% There are horrible people who, instead of solving a problem, tangle it up
% and make it harder to solve for anyone who wants to deal with it. Whoever
% does not know how to hit the nail on the head should be asked not to hit
% it at all.
%
% Friedrich Nietzsche

