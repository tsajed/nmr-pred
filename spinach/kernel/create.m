% Spin system and interaction specification.
%
% <http://spindynamics.org/wiki/index.php?title=Spin_system_specification>

function spin_system=create(sys,inter)

% Close all open files
fclose('all');

% Locate the root and run sanity checks
if isempty(which('existentials'))
    
    % Tell the user to RTFM
    error(['paths have not been correctly set - please follow installation '...
           'instructions (make sure you had included the subdirectories).']);
       
else
    
    % Set the root directory
    root_dir=which('existentials');
    spin_system.sys.root_dir=root_dir(1:end-32);
    
    % Run existential checks
    existentials();
    
end

% Validate input
grumble(sys,inter);

% Decide output destination
if isfield(sys,'output')&&strcmp(sys.output,'hush')
    
    % Hush the output
    spin_system.sys.output='hush';
    
elseif (~isfield(sys,'output'))||strcmp(sys.output,'console')
    
    % Print to the console
    spin_system.sys.output=1;
    
else
    
    % Print to a user-specified file
    spin_system.sys.output=fopen(sys.output,'a');
    
end

% Decide scratch destination
if isfield(sys,'scratch')
    
    % Scratch to a user-specified directory
    spin_system.sys.scratch=sys.scratch;
    
else
    
    % Scratch to the default directory
    spin_system.sys.scratch=[spin_system.sys.root_dir '/scratch'];
    
end

% Show version banner
banner(spin_system,'version_banner');

% Internal algorithms to disable
if isfield(sys,'disable')
    
    % Absorb user-specified values
    spin_system.sys.disable=sys.disable;
    
else
    
    % Disable nothing by default
    spin_system.sys.disable={};
    
end

% Internal algorithms to enable
if isfield(sys,'enable')
    
    % Absorb user-specified values
    spin_system.sys.enable=sys.enable;
    
else
    
    % Enable nothing by default
    spin_system.sys.enable={};
    
end

% Cut-offs and tolerances
spin_system=tolerances(spin_system,sys);

% Disabled features report
if ~isempty(spin_system.sys.disable)
    report(spin_system,'WARNING: the following functionality is disabled');
    if ismember('pt',spin_system.sys.disable),        report(spin_system,'         > automatic detection of non-interacting subspaces'); end
    if ismember('symmetry',spin_system.sys.disable),  report(spin_system,'         > permutation symmetry factorization'); end
    if ismember('krylov',spin_system.sys.disable),    report(spin_system,'         > Krylov propagation inside evolution() function'); end
    if ismember('clean-up',spin_system.sys.disable),  report(spin_system,'         > sparse array clean-up'); end
    if ismember('dss',spin_system.sys.disable),       report(spin_system,'         > destination state screening inside evolution() function'); end
    if ismember('expv',spin_system.sys.disable),      report(spin_system,'         > Krylov propagation inside step() function'); end
    if ismember('trajlevel',spin_system.sys.disable), report(spin_system,'         > trajectory analysis inside evolution() function'); end
    if ismember('merge',spin_system.sys.disable),     report(spin_system,'         > small subspace merging in the evolution() function'); end
    if ismember('norm_coil',spin_system.sys.disable), report(spin_system,'         > coil normalization in the evolution() function'); end
    if ismember('colorbar',spin_system.sys.disable),  report(spin_system,'         > colorbar drawing by plotting utilities'); end
end

% Enabled features report
if ~isempty(spin_system.sys.enable)
    report(spin_system,'WARNING: the following functionality is enabled');
    if ismember('gpu',spin_system.sys.enable),        report(spin_system,'         > GPU arithmetic'); end
    if ismember('caching',spin_system.sys.enable),    report(spin_system,'         > propagator caching'); end
    if ismember('greedy',spin_system.sys.enable),     report(spin_system,'         > greedy parallelisation'); end
    if ismember('xmemlist',spin_system.sys.enable),   report(spin_system,'         > state-cluster cross-membership list generation'); end
    if ismember('paranoia',spin_system.sys.enable),   report(spin_system,'         > paranoid accuracy settings'); end
end

% Control greedy parallelisation
warning('off','MATLAB:maxNumCompThreads:Deprecated');
if ismember('greedy',spin_system.sys.enable)
    spmd
        warning('off','MATLAB:maxNumCompThreads:Deprecated');
        maxNumCompThreads(feature('numcores'));
    end
    maxNumCompThreads(feature('numcores'));
else
    spmd
        warning('off','MATLAB:maxNumCompThreads:Deprecated');
        maxNumCompThreads(1);
    end
    maxNumCompThreads(feature('numcores'));
end

% Switch opengl to software on PCs
if ispc, opengl('software'); end

% GPU devices
if ismember('gpu',spin_system.sys.enable)
    
    % Inform the user
    report(spin_system,'looking for supported GPU devices...');
    
    % See what user says
    if isfield(sys,'gpuids')
        
        % Absorb user-specified values
        spin_system.sys.gpuids=sys.gpuids;
        
    else
        
        % Query the system
        spin_system.sys.gpuids=1:gpuDeviceCount;
        
    end
    
    % Clean up and inform the user
    if isempty(spin_system.sys.gpuids)
        spin_system.sys.enable=setdiff(spin_system.sys.enable,{'gpu'});
        report(spin_system,'WARNING - no CUDA capable GPUs detected');
    else  
        report(spin_system,['GPU devices to be used ' num2str(spin_system.sys.gpuids)]);
    end
    
else
    
    % Set GPU devices to none
    spin_system.sys.gpuids=[];
   
end

% Spin system banner
banner(spin_system,'spin_system_banner');

% Number and types of spins
spin_system.comp.isotopes=sys.isotopes;
spin_system.comp.nspins=numel(spin_system.comp.isotopes);

% Text labels for spins
if isfield(sys,'labels')
    
    % Get labels from the user
    spin_system.comp.labels=sys.labels;
    
else
    
    % Set labels to empty
    spin_system.comp.labels=cell(spin_system.comp.nspins,1);
    
end

% Multiplicities and magnetogyric ratios
spin_system.comp.mults=zeros(1,spin_system.comp.nspins);
spin_system.inter.gammas=zeros(1,spin_system.comp.nspins);
for n=1:spin_system.comp.nspins
    [spin_system.inter.gammas(n),spin_system.comp.mults(n)]=spin(sys.isotopes{n});
end
report(spin_system,['a total of ' num2str(spin_system.comp.nspins) ' particles in the simulation, of which '...
                                  num2str(nnz(spin_system.comp.mults==1)) ' have a zero spin.']);

% Order matrix
if isfield(inter,'order_matrix')
    spin_system.inter.order_matrix=inter.order_matrix;
else
    spin_system.inter.order_matrix=zeros(3);
end

% Primary magnet
spin_system.inter.magnet=sys.magnet;
report(spin_system,['magnetic induction of ' num2str(spin_system.inter.magnet,'%0.5g') ' Tesla ('...
                    num2str(-1e-6*spin_system.inter.magnet*spin('1H')/(2*pi),'%0.5g') ' MHz proton frequency, '...
                    num2str(-1e-9*spin_system.inter.magnet*spin('E' )/(2*pi),'%0.5g') ' GHz electron frequency).']);

% Compute carrier frequencies
spin_system.inter.basefrqs=-spin_system.inter.gammas*spin_system.inter.magnet;

% Preallocate Zeeman tensor array
spin_system.inter.zeeman.matrix=mat2cell(zeros(3*spin_system.comp.nspins,3),3*ones(spin_system.comp.nspins,1));

% Process Zeeman interactions
if isfield(inter,'zeeman')

    % Absorb eigenvalues and Euler angles
    if isfield(inter.zeeman,'eigs')
        for n=1:spin_system.comp.nspins
            if significant(inter.zeeman.eigs{n},0)
                if (~isfield(inter.zeeman,'euler'))||isempty(inter.zeeman.euler{n})
                    S=eye(3,3);
                else
                    S=euler2dcm(inter.zeeman.euler{n});
                end
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+S*diag(inter.zeeman.eigs{n})*transpose(S);
            end
        end
    end
    
    % Absorb tensors
    if isfield(inter.zeeman,'matrix')
        for n=1:spin_system.comp.nspins
            if significant(inter.zeeman.matrix{n},0)
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+inter.zeeman.matrix{n};
            end
        end
    end
    
    % Absorb scalars
    if isfield(inter.zeeman,'scalar')
         for n=1:spin_system.comp.nspins
            if significant(inter.zeeman.scalar{n},0)
                spin_system.inter.zeeman.matrix{n}=spin_system.inter.zeeman.matrix{n}+eye(3)*inter.zeeman.scalar{n};
            end
        end
    end

    % Report back to the user
    summary(spin_system,'zeeman','summary of Zeeman interactions as supplied (ppm for nuclei, g-tensor for electrons)');

    % Convert to angular frequencies
    for n=1:spin_system.comp.nspins
        switch spin_system.comp.isotopes{n}(1)
            case 'E'
                % For electrons, assume that the g-factor is given and compute the offset from the free electron g-factor
                spin_system.inter.zeeman.matrix{n}=(spin_system.inter.zeeman.matrix{n}-eye(3)*spin_system.tols.freeg)*(spin_system.inter.basefrqs(n)/spin_system.tols.freeg);
            otherwise
                % For nuclei, assume that the chemical shift is given and compute the corresponding offset 
                spin_system.inter.zeeman.matrix{n}=(1e-6)*spin_system.inter.zeeman.matrix{n}*spin_system.inter.basefrqs(n);
        end
    end
    
    % Clean up the result
    for n=1:spin_system.comp.nspins
        if negligible(spin_system.inter.zeeman.matrix{n},spin_system.tols.inter_cutoff)||(spin(spin_system.comp.isotopes{n})==0)
            spin_system.inter.zeeman.matrix{n}=[];
        end
    end
   
else
    
    % Warn the user that Zeeman tensors have not been found
    report(spin_system,'WARNING - no Zeeman interactions supplied, magnet frequencies assumed.');
    
end

% Absorb reaction rates
if isfield(inter,'chem')&&...
   isfield(inter.chem,'parts')&&...
   isfield(inter.chem,'rates')&&...
   isfield(inter.chem,'concs')
    
    % Assign the data structure
    spin_system.chem.parts=inter.chem.parts;
    spin_system.chem.rates=inter.chem.rates;
    spin_system.chem.concs=inter.chem.concs;
    
else
    
    % No chemical reactions and unit concentration
    spin_system.chem.parts={(1:spin_system.comp.nspins)};
    spin_system.chem.rates=0;
    spin_system.chem.concs=1;
    
end

% Absorb magnetization flux rates
if isfield(inter,'chem')&&...
   isfield(inter.chem,'flux_rate')&&...
   isfield(inter.chem,'flux_type')
    
    % Assign the data structure
    spin_system.chem.flux_rate=inter.chem.flux_rate;
    spin_system.chem.flux_type=inter.chem.flux_type;
    
else
    
    % No magnetization fluxes
    spin_system.chem.flux_rate=spalloc(spin_system.comp.nspins,spin_system.comp.nspins,0);
    spin_system.chem.flux_type='intramolecular';
    
end

% Report back to the user
summary(spin_system,'chemistry','chemical process summary');
    
% Preallocate coupling tensor array
report(spin_system,'initializing interaction arrays...');
spin_system.inter.coupling.matrix=mat2cell(zeros(3*spin_system.comp.nspins),3*ones(spin_system.comp.nspins,1),3*ones(spin_system.comp.nspins,1));

% Process coordinates
if isfield(inter,'coordinates')

    % Absorb coordinates into data structure
    spin_system.inter.coordinates=inter.coordinates;
    
    % Report back to the user
    summary(spin_system,'coordinates','atomic coordinates (Angstrom)');
    
    % Process periodic boundary conditions
    if isfield(inter,'pbc')
        
        % Absorb translation vectors into the data structure
        spin_system.inter.pbc=inter.pbc;
        
        % Report back to the user
        summary(spin_system,'pbc','PBC translation vectors (Angstrom)');
        
    else
        
        % Write an empty array
        spin_system.inter.pbc={};
        
        % Report back to the user
        report(spin_system,'periodic boundary conditions not specified, assuming a standalone spin system.');
        
    end
    
    % Call dipolar coupling module
    spin_system=dipolar(spin_system);
    
else
    
    % Warn the user that coordinates have not been found
    report(spin_system,'WARNING - no coordinates given, point dipolar interactions assumed to be zero.');
    
    % Set an empty coordinate array
    spin_system.inter.coordinates=cell(spin_system.comp.nspins,1);
    
    % Set proximity matrix to isolated spins
    spin_system.inter.proxmatrix=speye(spin_system.comp.nspins,spin_system.comp.nspins);
    
end

% Absorb user-specified couplings
if isfield(inter,'coupling')
    
    % Inform the user
    report(spin_system,'processing coupling data...');
    
    % Absorb eigenvalues and Euler angles
    if isfield(inter.coupling,'eigs')
        [rows,cols,~]=find(cellfun(@norm,inter.coupling.eigs)>spin_system.tols.inter_cutoff);
        for n=1:numel(rows)
            if isempty(inter.coupling.euler{rows(n),cols(n)})
                S=eye(3,3);
            else
                S=euler2dcm(inter.coupling.euler{rows(n),cols(n)});
            end
            spin_system.inter.coupling.matrix{rows(n),cols(n)}=...
            spin_system.inter.coupling.matrix{rows(n),cols(n)}+2*pi*S*diag(inter.coupling.eigs{rows(n),cols(n)})*S';
        end
    end
 
    % Absorb coupling tensors
    if isfield(inter.coupling,'matrix')
        [rows,cols,~]=find(cellfun(@norm,inter.coupling.matrix)>spin_system.tols.inter_cutoff);
        for n=1:numel(rows)
            spin_system.inter.coupling.matrix{rows(n),cols(n)}=...
            spin_system.inter.coupling.matrix{rows(n),cols(n)}+2*pi*inter.coupling.matrix{rows(n),cols(n)};
        end
    end
        
    % Absorb scalar couplings
    if isfield(inter.coupling,'scalar')
        [rows,cols,~]=find(cellfun(@norm,inter.coupling.scalar)>spin_system.tols.inter_cutoff);
        for n=1:numel(rows)
            spin_system.inter.coupling.matrix{rows(n),cols(n)}=...
            spin_system.inter.coupling.matrix{rows(n),cols(n)}+2*pi*eye(3)*inter.coupling.scalar{rows(n),cols(n)};
        end
    end
    
else
    
    % Warn the user that couplings have not been found
    report(spin_system,'WARNING - no couplings given, zeros assumed.');

end

% Order up coupling tensors
[rows,cols,~]=find(cellfun(@norm,spin_system.inter.coupling.matrix)>0);
for n=1:numel(rows)
    if rows(n)>cols(n)
        spin_system.inter.coupling.matrix{cols(n),rows(n)}=...
        spin_system.inter.coupling.matrix{cols(n),rows(n)}+spin_system.inter.coupling.matrix{rows(n),cols(n)};
        spin_system.inter.coupling.matrix{rows(n),cols(n)}=zeros(3);
    end
end

% Clean up the result
for n=1:spin_system.comp.nspins
    for k=1:spin_system.comp.nspins
        if negligible(spin_system.inter.coupling.matrix{n,k},spin_system.tols.inter_cutoff)||...
          (spin(spin_system.comp.isotopes{n})==0)||(spin(spin_system.comp.isotopes{k})==0)
            spin_system.inter.coupling.matrix{n,k}=[];
        end
    end
end

% Check for inter-subsystem couplings
for n=1:numel(spin_system.chem.parts)
    for k=1:numel(spin_system.chem.parts)
        if (n~=k)
            coupling_block=spin_system.inter.coupling.matrix(spin_system.chem.parts{n},spin_system.chem.parts{k});
            if ~all(cellfun(@isempty,coupling_block(:)))
                error('couplings detected between spins in different chemical species.');
            end
        end
    end
end

% Report back to the user
summary(spin_system,'couplings','summary of spin-spin couplings (angular frequencies)');

% Temperature
if ~isfield(inter,'temperature')
    
    % High-temperature approximation
    inter.temperature=[];
    
    % Print a warning to the user
    report(spin_system,'WARNING - high-temperature approximation');
    
else
    
    % Absorb the temperature specified
    spin_system.rlx.temperature=inter.temperature;
    
    % Report back to the user
    report(spin_system,['spin temperature: ' num2str(spin_system.rlx.temperature) ' Kelvin']);
    
end

% Relaxation superoperator
if isfield(inter,'relaxation')
    spin_system.rlx.theories=inter.relaxation;
    for n=1:numel(spin_system.rlx.theories)
        report(spin_system,['relaxation theory ' num2str(n) ': ' spin_system.rlx.theories{n}]);
    end
else
    spin_system.rlx.theories={};
    report(spin_system,'WARNING - no relaxation theory specified');
end

% Rotational correlation times
if isfield(inter,'tau_c')
    spin_system.rlx.tau_c=inter.tau_c;
    report(spin_system,['rotational correlation time(s): ' num2str(spin_system.rlx.tau_c) ' seconds']);
else
    spin_system.rlx.tau_c=0;
end

% Terms to keep in the relaxation superoperator
if ~isempty(spin_system.rlx.theories)
    spin_system.rlx.keep=inter.rlx_keep;
    report(spin_system,['terms to keep in the relaxation superoperator: ' spin_system.rlx.keep]);
end

% The fate of the dynamic frequency shift
if isfield(inter,'rlx_dfs')
    spin_system.rlx.dfs=inter.rlx_dfs;
else
    spin_system.rlx.dfs='ignore';
end
if ~isempty(spin_system.rlx.theories)
    report(spin_system,['action to take on dynamic frequency shifts: ' spin_system.rlx.dfs]);
end

% SRFK correlation time
if isfield(inter,'srfk_tau_c')
    spin_system.rlx.srfk_tau_c=inter.srfk_tau_c;
else
    spin_system.rlx.srfk_tau_c=0;
end
if ismember('SRFK',spin_system.rlx.theories)
    report(spin_system,['SRFK correlation time: ' num2str(spin_system.rlx.srfk_tau_c) ' seconds']);
end

% SRFK assumptions
if isfield(inter,'srfk_assume')
    spin_system.rlx.srfk_assume=inter.srfk_assume;
else
    spin_system.rlx.srfk_assume='';
end
if ismember('SRFK',spin_system.rlx.theories)
    report(spin_system,['SRFK Hamiltonian assumptions: ' spin_system.rlx.srfk_assume]);
end

% SRFK modulation depths
if isfield(inter,'srfk_mdepth')
    spin_system.rlx.srfk_mdepth=inter.srfk_mdepth;
else
    spin_system.rlx.srfk_mdepth=[];
end
if ismember('SRFK',spin_system.rlx.theories)&&(~isempty(spin_system.rlx.srfk_mdepth))
    for n=1:spin_system.comp.nspins
        for k=1:spin_system.comp.nspins
            if ~isempty(spin_system.rlx.srfk_mdepth{n,k})
                report(spin_system,['SRFK modulation depth for spins ' ...
                                    num2str(n) ',' num2str(k) ': ' ...
                                    num2str(spin_system.rlx.srfk_mdepth{n,k}) ' Hz']);
            end
        end
    end
end

% Equilibrium state
if ~isempty(spin_system.rlx.theories)
    
    % Absorb user setting
    spin_system.rlx.equilibrium=inter.equilibrium;
    
    % Print back the notice
    report(spin_system,['thermalisation method: ' spin_system.rlx.equilibrium]);
    
end

% Isotropic damping
if isfield(inter,'damp_rate')
    spin_system.rlx.damp_rate=inter.damp_rate;
    report(spin_system,['isotropic damping rate: ' num2str(spin_system.rlx.damp_rate) ' Hz.']);
else
    spin_system.rlx.damp_rate=[];
end

% Anisotropic damping
if isfield(inter,'hstrain_rates')
    spin_system.rlx.hstrain_rates=inter.hstrain_rates;
    report(spin_system,['anisotropic damping rates: ' num2str(spin_system.rlx.hstrain_rates) ' Hz.']);
else
    spin_system.rlx.hstrain_rates=[];
end

% User-supplied R1 relaxation rates for T1/T2 model
if isfield(inter,'r1_rates')
    spin_system.rlx.r1_rates=inter.r1_rates;
else
    spin_system.rlx.r1_rates=[];
end

% User-supplied R2 relaxation rates for T1/T2 model
if isfield(inter,'r2_rates')
    spin_system.rlx.r2_rates=inter.r2_rates;
else
    spin_system.rlx.r2_rates=[];
end

% User-supplied R1 relaxation rates for Lindblad theory
if isfield(inter,'lind_r1_rates')
    spin_system.rlx.lind_r1_rates=inter.lind_r1_rates;
else
    spin_system.rlx.lind_r1_rates=[];
end

% User-supplied R2 relaxation rates for Lindblad theory
if isfield(inter,'lind_r2_rates')
    spin_system.rlx.lind_r2_rates=inter.lind_r2_rates;
else
    spin_system.rlx.lind_r2_rates=[];
end

% User-supplied Weizmann R1e relaxation rates
if isfield(inter,'weiz_r1e')
    spin_system.rlx.weiz_r1e=inter.weiz_r1e;
else
    spin_system.rlx.weiz_r1e=[];
end

% User-supplied Nottingham R1e relaxation rates
if isfield(inter,'nott_r1e')
    spin_system.rlx.nott_r1e=inter.nott_r1e;
else
    spin_system.rlx.nott_r1e=[];
end

% User-supplied Weizmann R1n relaxation rates
if isfield(inter,'weiz_r1n')
    spin_system.rlx.weiz_r1n=inter.weiz_r1n;
else
    spin_system.rlx.weiz_r1n=[];
end

% User-supplied Nottingham R1n relaxation rates
if isfield(inter,'nott_r1n')
    spin_system.rlx.nott_r1n=inter.nott_r1n;
else
    spin_system.rlx.nott_r1n=[];
end

% User-supplied Weizmann R1d relaxation rates
if isfield(inter,'weiz_r1d')
    spin_system.rlx.weiz_r1d=inter.weiz_r1d;
else
    spin_system.rlx.weiz_r1d=[];
end

% User-supplied Weizmann R2e relaxation rates
if isfield(inter,'weiz_r2e')
    spin_system.rlx.weiz_r2e=inter.weiz_r2e;
else
    spin_system.rlx.weiz_r2e=[];
end

% User-supplied Nottingham R2e relaxation rates
if isfield(inter,'nott_r2e')
    spin_system.rlx.nott_r2e=inter.nott_r2e;
else
    spin_system.rlx.nott_r2e=[];
end

% User-supplied Weizmann R2n relaxation rates
if isfield(inter,'weiz_r2n')
    spin_system.rlx.weiz_r2n=inter.weiz_r2n;
else
    spin_system.rlx.weiz_r2n=[];
end

% User-supplied Nottingham R2n relaxation rates
if isfield(inter,'nott_r2n')
    spin_system.rlx.nott_r2n=inter.nott_r2n;
else
    spin_system.rlx.nott_r2n=[];
end

% User-supplied Weizmann R2d relaxation rates
if isfield(inter,'weiz_r2d')
    spin_system.rlx.weiz_r2d=inter.weiz_r2d;
else
    spin_system.rlx.weiz_r2d=[];
end

% Report relaxation rates back to the user
if isfield(inter,'r1_rates')&&isfield(inter,'r2_rates')
    summary(spin_system,'rlx_rates_t1_t2','relaxation rates (Hz) for T1/T2 theory');
end
if isfield(inter,'lind_r1_rates')&&isfield(inter,'lind_r2_rates')
    summary(spin_system,'rlx_rates_lindblad','relaxation rates (Hz) for Lindblad theory');
end
if isfield(inter,'nott_r1e')&&isfield(inter,'nott_r2e')&&...
   isfield(inter,'nott_r1n')&&isfield(inter,'nott_r2n')
    summary(spin_system,'rlx_rates_nott','relaxation rates (Hz) for Nottingham DNP theory');
end
if isfield(inter,'weiz_r1e')&&isfield(inter,'weiz_r2e')&&...
   isfield(inter,'weiz_r1n')&&isfield(inter,'weiz_r2n')&&...
   isfield(inter,'weiz_r1d')&&isfield(inter,'weiz_r2d')
    summary(spin_system,'rlx_rates_weiz','relaxation rates (Hz) for Weizmann DNP theory');
end

% Absorb radical recombination parameters
if isfield(inter,'chem')&&isfield(inter.chem,'rp_theory')

    % Absorb theory
    spin_system.chem.rp_theory=inter.chem.rp_theory;
    report(spin_system,['radical recombination theory set to ' spin_system.chem.rp_theory]);

    % Absorb spins
    spin_system.chem.rp_electrons=inter.chem.rp_electrons;
    report(spin_system,['recombining electrons at positions ' num2str(spin_system.chem.rp_electrons)]);

    % Absorb rates
    spin_system.chem.rp_rates=inter.chem.rp_rates;
    report(spin_system,['singlet recombination rate ' num2str(spin_system.chem.rp_rates(1)) ' Hz.']);
    report(spin_system,['triplet recombination rate ' num2str(spin_system.chem.rp_rates(2)) ' Hz.']);

else
    
    spin_system.chem.rp_theory='';
    spin_system.chem.rp_electrons=[];
    spin_system.chem.rp_rates=[];
    
end

end

% Input validation function
function grumble(sys,inter)

% Check Matlab component versions
if verLessThan('matlab','9.0')
    error('Spinach requires 64-bit Matlab R2016a or later.');
end
if verLessThan('distcomp','6.8')
    error('Spinach requires Matlab Parallel Computing Toolbox version 6.8 or later.');
end
if verLessThan('optim','7.4')
    error('Spinach requires Matlab Optimization Toolbox version 7.4 or later.');
end

% Check the output switch
if isfield(sys,'output')
    if ~ischar(sys.output)
        error('sys.output must be a character string.');
    end
end

% Check the scratch folder
if isfield(sys,'scratch')
    if ~ischar(sys.scratch)
        error('sys.scratch must be a character string.');
    end
    if ~exist(sys.scratch,'dir')
        error('the specified scratch directory does not exist.');
    end
end

% Check the disable switch
if isfield(sys,'disable')
    if (~iscell(sys.disable))||any(~cellfun(@ischar,sys.disable))
        error('sys.disable must be a cell array of strings.');
    end
    if any(~ismember(sys.disable,{'zte','pt','symmetry','krylov','clean-up','dss','expv','trajlevel','merge','norm_coil','colorbar'}))
        error(['allowed values for sys.disable are ''zte'', ''pt'', ''symmetry'', ''krylov'', ''clean-up'','...
               ' ''dss'', ''expv'', ''merge'', ''colorbar'', ''norm_coil'' and ''trajlevel''.']);
    end
end

% Check the enable switch
if isfield(sys,'enable')
    if (~iscell(sys.enable))||any(~cellfun(@ischar,sys.enable))
        error('sys.enable must be a cell array of strings.');
    end
    if any(~ismember(sys.enable,{'gpu','caching','xmemlist','greedy','paranoia','cowboy'}))
        error('allowed values for sys.enable are ''gpu'', ''xmemlist'', ''greedy'', ''paranoia'', ''cowboy'' and ''caching''.');
    end
end

% Check GPU parameters
if isfield(sys,'gpuids')
    
    % Check the specification
    if ~isnumeric(sys.gpuids)||any(mod(sys.gpuids,1)~=0)
        error('sys.gpuids must be a vector of integers.');
    end
    
    % Check if the devices are there
    if any(sys.gpuids>gpuDeviceCount)
        error('GPU device with the specified ID does not exist.');
    end
    
end

% Check isotopes variable
if ~isfield(sys,'isotopes')
    error('sys.isotopes subfield must be present.');
elseif isempty(sys.isotopes)
    error('sys.isotopes cell array cannot be empty.');
elseif ~iscell(sys.isotopes)
    error('sys.isotopes must be a cell array.');
elseif ~all(cellfun(@ischar,sys.isotopes))
    error('all elements of sys.isotopes cell array must be character strings.');
end

% Check labels variable
if isfield(sys,'labels')
    if ~iscell(sys.labels)
        error('sys.labels must be a cell array.');
    elseif isempty(sys.labels)
        error('sys.labels cell array cannot be empty.');
    elseif ~all(cellfun(@ischar,sys.labels))
        error('all elements of sys.labels cell array must be character strings.');
    elseif numel(sys.labels)~=numel(sys.isotopes)
        error('the length of sys.labels must be the same as the length of sys.isotopes.');
    end
end

% Check the order matrix
if isfield(inter,'order_matrix')
    
    % Make sure we have a 3x3 matrix
    if (~isnumeric(inter.order_matrix))||(~ismatrix(inter.order_matrix))||...
       (any(size(inter.order_matrix)~=[3 3]))||(~isreal(inter.order_matrix))
        error('inter.order_matrix must be a 3x3 matrix of real numbers.');
    end
    
    % Make sure the matrix is traceless
    if trace(inter.order_matrix)>10*eps()
        error('inter.order_matrix must be traceless.');
    end
    
    % Make sure the matrix is symmetric
    if norm(inter.order_matrix-transpose(inter.order_matrix))>10*eps()
        error('inter.order_matrix must be symmetric.');
    end
    
end

% Check magnet induction
if ~isfield(sys,'magnet')
    error('magnet induction must be specified in sys.magnet varaible.');
else
    if (~isnumeric(sys.magnet))||(~isreal(sys.magnet))||(numel(sys.magnet)~=1)
        error('sys.magnet must be a real number.');
    end
end

% Check Zeeman interactions
if isfield(inter,'zeeman')
    
    % Check eigenvalues / Euler angles specification
    if isfield(inter.zeeman,'eigs')
        
        % Check type
        if (~iscell(inter.zeeman.eigs))||any(~cellfun(@isnumeric,inter.zeeman.eigs))
            error('inter.zeeman.eigs must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Check dimensions
        if numel(inter.zeeman.eigs)~=numel(sys.isotopes)
            error('the number of elements in inter.zeeman.eigs must match the number of spins.');
        end
        
        % Make sure eulers exist
        if ~isfield(inter.zeeman,'euler')
            error('inter.zeeman.euler variable must be set together with inter.zeeman.eigs.');
        end
        
        % Make sure eulers are cells
        if (~iscell(inter.zeeman.euler))||any(~cellfun(@isnumeric,inter.zeeman.euler))
            error('inter.zeeman.euler must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Make sure eulers have the correct length
        if ~all(size(inter.zeeman.eigs)==size(inter.zeeman.euler))
            error('inter.zeeman.eigs and inter.zeeman.euler variables must have the same dimension.');
        end
        
        % Make sure all non-empty elements are real 1x3 vectors
        for n=1:numel(sys.isotopes)
            
            % For eigenvalues
            if ~isempty(inter.zeeman.eigs{n})
                if (~all(size(inter.zeeman.eigs{n})==[1 3]))||...
                   (~isnumeric(inter.zeeman.eigs{n}))||...
                   (~isreal(inter.zeeman.eigs{n}))
                    error('non-empty elements of inter.zeeman.eigs must be real 1x3 vectors.');
                end
            end
            
            % For Euler angles
            if ~isempty(inter.zeeman.euler{n})
                if (~all(size(inter.zeeman.euler{n})==[1 3]))||...
                   (~isnumeric(inter.zeeman.euler{n}))||...
                   (~isreal(inter.zeeman.euler{n}))
                    error('non-empty elements of inter.zeeman.euler must be real 1x3 vectors.');
                end
            end
            
            % For simultaneity
            if (isempty(inter.zeeman.eigs{n})&&(~isempty(inter.zeeman.euler{n})))||...
               (isempty(inter.zeeman.euler{n})&&(~isempty(inter.zeeman.eigs{n})))
                error('inter.zeeman.eigs and inter.zeeman.euler must have identical non-empty cell patterns.');
            end
            
        end
        
    end
    
    % Check matrix specification
    if isfield(inter.zeeman,'matrix')
        
        % Check type
        if ~iscell(inter.zeeman.matrix)
            error('inter.zeeman.matrix must be a cell array of empty matrices or 3x3 matrices.');
        elseif size(inter.zeeman.matrix,1)~=1
            error('inter.zeeman.matrix cell array must have dimension 1 x nspins.');
        elseif ~all(cellfun(@isnumeric,inter.zeeman.matrix))
            error('all elements of inter.zeeman.matrix cell array must be numeric.');
        end 
        
        % Check length
        if numel(inter.zeeman.matrix)~=numel(sys.isotopes)
            error('the number of elements in the inter.zeeman.matrix array should match the number of spins.');
        end
        
        % Make sure all non-empty elements are real 3x3 matrices
        for n=1:numel(sys.isotopes)
            if ~isempty(inter.zeeman.matrix{n})
                if (~all(size(inter.zeeman.matrix{n})==[3 3]))||...
                   (~isnumeric(inter.zeeman.matrix{n}))||...
                   (~isreal(inter.zeeman.matrix{n}))
                    error('non-empty elements of inter.zeeman.matrix must be real 3x3 matrices.');
                end
            end
        end
        
    end
    
    % Check scalars
    if isfield(inter.zeeman,'scalar')
        
        % Check type
        if ~iscell(inter.zeeman.scalar)||any(~cellfun(@isnumeric,inter.zeeman.scalar))
            error('inter.zeeman.scalar must be a cell array of empty matrices or 1x1 matrices.');
        end
        
        % Check length
        if numel(inter.zeeman.scalar)~=numel(sys.isotopes)
            error('the number of elements in the inter.zeeman.scalar array should match the number of spins.');
        end
        
        % Make sure all non-empty elements are real numbers
        for n=1:numel(sys.isotopes)
            if ~isempty(inter.zeeman.scalar{n})
                if (numel(inter.zeeman.scalar{n})~=1)||...
                   (~isnumeric(inter.zeeman.scalar{n}))||...
                   (~isreal(inter.zeeman.scalar{n}))
                    error('non-empty elements of inter.zeeman.scalar must be numbers.');
                end
            end
        end
        
    end
    
end

% Check coordinates
if isfield(inter,'coordinates')
    
    % Check type
    if ~iscell(inter.coordinates)||any(~cellfun(@isnumeric,inter.coordinates))
        error('inter.coordinates must be a cell array of empty vectors or 1x3 vectors.');
    end
    
    % Check size
    if numel(inter.coordinates)~=numel(sys.isotopes)
        error('the number of elements in inter.coordinates must match the number of spins.')
    end
    
    % Check contents
    for n=1:numel(sys.isotopes)
        
        % Make sure we have real 3-vectors
        if ~isempty(inter.coordinates{n})
            if (~all(size(inter.coordinates{n})==[1 3]))||...
               (~isnumeric(inter.coordinates{n}))||...
               (~isreal(inter.coordinates{n}))
                error('non-empty elements of inter.coordinates must be real 1x3 vectors.');
            end
        end
        
    end
    
end

% Check periodic boundary conditions
if isfield(inter,'pbc')
    
    % Check type
    if ~iscell(inter.pbc)||any(~cellfun(@isvector,inter.pbc))
        error('inter.pbc must be a cell array of row vectors.');
    end
    
    % Check element numbers
    if ~ismember(numel(inter.pbc),[0 1 2 3])
        error('inter.pbc cell array must have zero, one, two or three row vectors as elements.');
    end
    
    % Check vector dimensions
    for n=1:numel(inter.pbc)
        if (~all(size(inter.pbc{n})==[1 3]))||(~isreal(inter.pbc{n}))
            error('all elements of inter.pbc.cell array must be row vectors with three real elements.');
        end
    end
    
    % Check linear independence
    if ((numel(inter.pbc)==2)&&rank([inter.pbc{1}; inter.pbc{2}],1e-3)<2)||...
       ((numel(inter.pbc)==3)&&rank([inter.pbc{1}; inter.pbc{2}; inter.pbc{3}],1e-3)<3)
        error('the vectors supplied in inter.pbc must be lineary independent.');
    end
    
end
            
% Check couplings
if isfield(inter,'coupling')
    
    % Check eigenvalues / Euler angles specification
    if isfield(inter.coupling,'eigs')
        
        % Check type
        if (~iscell(inter.coupling.eigs))||any(any(~cellfun(@isnumeric,inter.coupling.eigs)))
            error('inter.coupling.eigs must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Check dimensions
        if ~all(size(inter.coupling.eigs)==[numel(sys.isotopes) numel(sys.isotopes)])
            error('both dimensions of inter.coupling.eigs must match the number of spins.');
        end
        
        % Check quadratic couplings
        for n=1:numel(sys.isotopes)
            [~,mult]=spin(sys.isotopes{n});
            if (norm(inter.coupling.eigs{n,n})>0)&&(mult<3)
                error('quadratic couplings cannot be specified for spin-1/2 particles.');
            elseif abs(sum(inter.coupling.eigs{n,n}))>1e-6
                error('quadratic couplings cannot have a non-zero trace.');
            end
        end
        
        % Make sure eulers exist
        if ~isfield(inter.coupling,'euler')
            error('inter.coupling.euler array must be set together with inter.coupling.eigs.');
        end
        
        % Make sure eulers are cells
        if ~iscell(inter.coupling.euler)||any(any(~cellfun(@isnumeric,inter.coupling.euler)))
            error('inter.coupling.euler must be a cell array of empty vectors or 1x3 vectors.');
        end
        
        % Make sure eulers have the correct length
        if ~all(size(inter.coupling.eigs)==size(inter.coupling.euler))
            error('inter.coupling.eigs and inter.coupling.euler arrays must have the same dimension.');
        end
        
        % Make sure all non-empty elements are real 1x3 vectors
        for n=1:numel(sys.isotopes)
            for k=1:numel(sys.isotopes)
            
                % For eigenvalues
                if ~isempty(inter.coupling.eigs{n,k})
                    if (~all(size(inter.coupling.eigs{n,k})==[1 3]))||...
                       (~isnumeric(inter.coupling.eigs{n,k}))||...
                       (~isreal(inter.coupling.eigs{n,k}))
                        error('non-empty elements of inter.coupling.eigs must be real 1x3 vectors.');
                    end
                end
                
                % For Euler angles
                if ~isempty(inter.coupling.euler{n,k})
                    if (~all(size(inter.coupling.euler{n,k})==[1 3]))||...
                       (~isnumeric(inter.coupling.euler{n,k}))||...
                       (~isreal(inter.coupling.euler{n,k}))
                        error('non-empty elements of inter.coupling.euler must be real 1x3 vectors.');
                    end
                end
                
                % For simultaneity
                if (isempty(inter.coupling.eigs{n,k})&&(~isempty(inter.coupling.euler{n,k})))||...
                   (isempty(inter.coupling.euler{n,k})&&(~isempty(inter.coupling.eigs{n,k})))
                    error('inter.coupling.eigs and inter.coupling.euler must have identical non-empty cell patterns.');
                end
                
            end
        end
        
    end
    
    % Check matrix specification
    if isfield(inter.coupling,'matrix')
        
        % Check type
        if ~iscell(inter.coupling.matrix)||any(any(~cellfun(@isnumeric,inter.coupling.matrix)))
            error('inter.coupling.matrix must be a cell array of empty matrices or 3x3 matrices.');
        end
        
        % Check dimensions
        if ~all(size(inter.coupling.matrix)==[numel(sys.isotopes) numel(sys.isotopes)])
            error('both dimensions of inter.coupling.matrix cell array should match the number of spins.');
        end
        
        % Check quadratic couplings
        for n=1:numel(sys.isotopes)
            [~,mult]=spin(sys.isotopes{n});
            if (norm(inter.coupling.matrix{n,n})>0)&&(mult<3)
                error('quadratic couplings cannot be specified for spin-1/2 particles.');
            elseif abs(trace(inter.coupling.matrix{n,n}))>1e-6
                error('quadratic couplings cannot have a non-zero trace.');
            end
        end
        
        % Make sure all non-empty elements are real 3x3 matrices
        for n=1:numel(sys.isotopes)
            for k=1:numel(sys.isotopes)
                if ~isempty(inter.coupling.matrix{n,k})
                    if (~all(size(inter.coupling.matrix{n,k})==[3 3]))||...
                       (~isnumeric(inter.coupling.matrix{n,k}))||...
                       (~isreal(inter.coupling.matrix{n,k}))
                        error('non-empty elements of inter.coupling.matrix must be real 3x3 matrices.');
                    end
                end
            end
        end
        
    end
    
    % Check scalars
    if isfield(inter.coupling,'scalar')
        
        % Check type
        if ~iscell(inter.coupling.scalar)||any(any(~cellfun(@isnumeric,inter.coupling.scalar)))
            error('inter.coupling.scalar must be a cell array of empty matrices or 1x1 matrices.');
        end
        
        % Check dimensions
        if ~all(size(inter.coupling.scalar)==[numel(sys.isotopes) numel(sys.isotopes)])
            error('both dimensions of inter.coupling.scalar array should match the number of spins.');
        end
        
        % Make sure all non-empty elements are real numbers
        if any(nonzeros(cellfun(@numel,inter.coupling.scalar))~=1)||...
           any(nonzeros(~cellfun(@isnumeric,inter.coupling.scalar)))||...
           any(nonzeros(~cellfun(@isreal,inter.coupling.scalar)))
            error('non-empty elements of inter.coupling.scalar must be real numbers.');
        end
        
        % Disallow quadratic scalar couplings
        for n=1:numel(sys.isotopes)
            if norm(inter.coupling.scalar{n,n})~=0
                error('scalar couplings cannot be quadratic.');
            end
        end
         
    end
    
end

% Check temperature - negative and complex values allowed
if isfield(inter,'temperature')&&(~isempty(inter.temperature))
    
    % Check type and dimension
    if (~isnumeric(inter.temperature))||(numel(inter.temperature)~=1)
        error('inter.tempearture must be a number.');
    end
    
end

% Check relaxation theory
if isfield(inter,'relaxation')
    
    % Check type
    if (~iscell(inter.relaxation))||any(~cellfun(@ischar,inter.relaxation))
        error('inter.relaxation must be a cell array of strings.');
    end
    
    % Check relaxation theories
    if ~all(ismember(inter.relaxation,{'damp','hstrain','t1_t2','redfield','lindblad',...
                                       'nottingham','weizmann','SRFK','SRSK'}))
        error('unrecognised relaxation theory specification.');
    end
    
    % Enforce term retention policy specification
    if ~isfield(inter,'rlx_keep')
        error('relaxation superoperator term retention policy must be specified in inter.rlx_keep field.');
    end
    
    % Enforce relaxation destination
    if ~isfield(inter,'equilibrium')
        error('relaxation destination must be specified in inter.equilibrium variable.');
    end
    
    % Enforce correlation time with Redfield theory
    if ismember('redfield',inter.relaxation)&&(~isfield(inter,'tau_c'))
        error('correlation time(s) must be specified with Redfield theory.');
    end
    
    % Enforce damping rate with isotropic damping
    if ismember('damp',inter.relaxation)&&(~isfield(inter,'damp_rate'))
        error('damping rate must be specified with non-selective damping.');
    end
    
    % Enforce damping rate with anisotropic damping
    if ismember('hstrain',inter.relaxation)&&(~isfield(inter,'hstrain_rates'))
        error('damping rates must be specified with anisotropic damping.');
    end
    
    % Enforce R1 and R2 rates with T1,T2 approximation
    if ismember('t1_t2',inter.relaxation)&&((~isfield(inter,'r1_rates'))||(~isfield(inter,'r2_rates')))
        error('R1 and R2 rates must be specified with extended T1,T2 relaxation theory.');
    end
    
    % Enforce R1 and R2 rates with Lindblad theory
    if ismember('lindblad',inter.relaxation)&&((~isfield(inter,'lind_r1_rates'))||(~isfield(inter,'lind_r2_rates')))
        error('R1 and R2 rates must be specified with Lindblad relaxation theory.');
    end
    
    % Enforce R1e, R2e, R1n and R2n rates with Nottingham DNP theory
    if ismember('nottingham',inter.relaxation)&&...
       ((~isfield(inter,'nott_r1e'))||(~isfield(inter,'nott_r2e'))||...
        (~isfield(inter,'nott_r1n'))||(~isfield(inter,'nott_r2n')))
        error(['R1e, R2e, R1n and R2n rates must be specified with '...
               'Nottingham DNP relaxation theory.']);
    end
    
    % Enforce R1e, R2e, R1n and R2n rates with Weizmann DNP theory
    if ismember('weizmann',inter.relaxation)&&...
       ((~isfield(inter,'weiz_r1e'))||(~isfield(inter,'weiz_r2e'))||...
        (~isfield(inter,'weiz_r1n'))||(~isfield(inter,'weiz_r2n')))
        error(['R1e, R2e, R1n and R2n rates must be specified with '...
               'Weizmann DNP relaxation theory.']);
    end
    
    % Enforce R1d and R2d with Weizmann DNP theory
    if ismember('weizmann',inter.relaxation)&&...
       ((~isfield(inter,'weiz_r1d'))||(~isfield(inter,'weiz_r2d')))
        error('R1d and R2d rates must be specified with Weizmann DNP relaxation theory.');
    end
    
    % Enforce two electrons with Nottingham DNP theory
    if ismember('nottingham',inter.relaxation)&&(nnz(strcmp('E',sys.isotopes))~=2)
        error('Nottingham DNP relaxation theory requires two electrons.');
    end
    
    % Enforce correlation time with SRFK
    if ismember('SRFK',inter.relaxation)&&(~isfield(inter,'srfk_tau_c'))
        error('SRFK requires modulation correlation time to be specified in inter.srfk_tau_c variable.');
    end
    
    % Enforce modulation depth with SRFK
    if ismember('SRFK',inter.relaxation)&&(~isfield(inter,'srfk_mdepth'))
        error('SRFK requires modulation depths to be specified in inter.srfk_mdepth variable.');
    end
    
    % Enforce specific assumptions with SRFK
    if ismember('SRFK',inter.relaxation)&&(~isfield(inter,'srfk_assume'))
        error('SRFK requires assumptions to be specified in inter.srfk_assume variable.');
    end
    
end

% Check rotational correlation time
if isfield(inter,'tau_c')
    
    % Check type and dimension
    if (~isnumeric(inter.tau_c))||(numel(inter.tau_c)>3)||(numel(inter.tau_c)==0)
        error('inter.tau_c must be a vector of size 1, 2 or 3.');
    end
    
    % Check value
    if (~isreal(inter.tau_c))||any(inter.tau_c<0)
        error('inter.tau_c must have non-negative real elements.');
    end
    
    % Enforce Redfield theory if tau_c is specified
    if (~isfield(inter,'relaxation'))||(~ismember('redfield',inter.relaxation))
        error('inter.tau_c requires Redfield relaxation theory.');
    end
    
end

% Check SRFK correlation time
if isfield(inter,'srfk_tau_c')
    
    % Check type and dimension
    if (~isnumeric(inter.srfk_tau_c))||(~isscalar(inter.srfk_tau_c))||(~isreal(inter.tau_c))||(inter.tau_c<0)
        error('inter.srfk_tau_c must be a non-negative real number.');
    end
    
    % Enforce SRFK theory if srfk_tau_c is specified
    if (~isfield(inter,'relaxation'))||(~ismember('SRFK',inter.relaxation))
        error('inter.srfk_tau_c requires SRFK relaxation theory.');
    end
    
end

% Check term retention
if isfield(inter,'rlx_keep')
    
    % Check type
    if ~ischar(inter.rlx_keep)
        error('inter.rlx_keep must be a string.');
    end
    
    % Check contents
    if ~ismember(inter.rlx_keep,{'diagonal','kite','secular','labframe'})
        error('allowed values for inter.rlx_keep are ''diagonal'', ''kite'', ''secular'' and ''labframe''.');
    end
    
end

% Check dynamic frequency shift retention
if isfield(inter,'rlx_dfs')
    
    % Check type
    if ~ischar(inter.rlx_dfs)
        error('inter.rlx_dfs must be a string.');
    end
    
    % Check contents
    if ~ismember(inter.rlx_dfs,{'keep','ignore'})
        error('allowed values for inter.rlx_dfs are ''keep'' and ''ignore''.');
    end
    
end
    
% Check SRFK assumptions
if isfield(inter,'srfk_assume')
    
    % Check type
    if ~ischar(inter.srfk_assume)
        error('inter.srfk_assume must be a string.');
    end
    
    % Enforce SRFK if SRFK assumptions are specified
    if (~isfield(inter,'relaxation'))||(~ismember('SRFK',inter.relaxation))
        error('inter.srfk_assume requires SRFK relaxation theory.');
    end
    
end

% Check SRFK modulation depths
if isfield(inter,'srfk_mdepth')

    % Check type
    if (~iscell(inter.srfk_mdepth))||(~all(cellfun(@isnumeric,inter.srfk_mdepth(:))))
        error('inter.srfk_mdepth must be a cell array of empty matrices or scalars.');
    end
    
    % Check dimensions
    if ~all(size(inter.srfk_mdepth)==[numel(sys.isotopes) numel(sys.isotopes)])
        error('both dimensions of inter.srfk_mdepth array should match the number of spins.');
    end
    
    % Make sure all non-empty elements are non-negative real numbers
    [rows,cols]=find(~cellfun(@isempty,inter.srfk_mdepth));
    for n=1:numel(rows)
        if (~isreal(inter.srfk_mdepth{rows(n),cols(n)}))||(inter.srfk_mdepth{rows(n),cols(n)}<0)
            error('non-empty elements of inter.srfk_mdepth must be non-negative real numbers.');
        end
    end
    
    % Disallow quadratic scalar couplings
    for n=1:numel(sys.isotopes)
        if norm(inter.srfk_mdepth{n,n})~=0
            error('scalar couplings cannot be quadratic.');
        end
    end
    
    % Enforce SRFK if SRFK modulation depths are specified
    if (~isfield(inter,'relaxation'))||(~ismember('SRFK',inter.relaxation))
        error('inter.srfk_mdepth requires SRFK relaxation theory.');
    end
    
end

% Check equilibrium switch
if isfield(inter,'equilibrium')
    
    % Check type
    if ~ischar(inter.equilibrium)
        error('inter.equilibrium must be a string.');
    end
    
    % Check contents
    if ~ismember(inter.equilibrium,{'zero','levante','dibari'})
        error('allowed values for inter.equilibrium are ''zero'', ''levante'' and ''dibari''.');
    end
    
end

% Check isotropic damping rate
if isfield(inter,'damp_rate')
    
    % Check value
    if (~isnumeric(inter.damp_rate))||(numel(inter.damp_rate)~=1)||...
       (~isreal(inter.damp_rate))||(inter.damp_rate<0)
        error('inter.damp_rate must be a non-negative real number.');
    end
    
    % Enforce isotropic damping if damp_rate is specified
    if ~ismember('damp',inter.relaxation)
        error('inter.damp_rate can only be specified with damp relaxation theory.');
    end
    
end

% Check anisotropic damping rates
if isfield(inter,'hstrain_rates')
    
    % Check values
    if (~isnumeric(inter.hstrain_rates))||(numel(inter.hstrain_rates)~=3)||...
       (~isreal(inter.hstrain_rates))||any(inter.hstrain_rates<0)
        error('inter.hstrain_rates must contain three non-negative real numbers.');
    end
    
    % Enforce anisotropic damping if damp_rate is specified
    if ~ismember('hstrain',inter.relaxation)
        error('inter.hstrain_rates can only be specified with damp relaxation theory.');
    end
    
end

% Check R1 rates for T1/T2
if isfield(inter,'r1_rates')
    
    % Check type and values
    if (~isvector(inter.r1_rates))||(~isnumeric(inter.r1_rates))||...
         any(inter.r1_rates<0)||any(~isreal(inter.r1_rates))
        error('inter.r1_rates must be a vector of positive real numbers.');
    end
    
    % Check dimension
    if numel(inter.r1_rates)~=numel(sys.isotopes)
        error('the number of elements in inter.r1_rates must be equal to the number of spins.');
    end
    
    % Enforce T1,T2 theory if inter.r1_rates rates are specified
    if ~ismember('t1_t2',inter.relaxation)
        error('inter.r1_rates can only be specified with T1,T2 relaxation theory.');
    end
    
end

% Check R1 rates for Lindblad
if isfield(inter,'lind_r1_rates')
    
    % Check type and values
    if (~isvector(inter.lind_r1_rates))||(~isnumeric(inter.lind_r1_rates))||...
         any(inter.lind_r1_rates<0)||any(~isreal(inter.lind_r1_rates))
        error('inter.lind_r1_rates must be a vector of positive real numbers.');
    end
    
    % Check dimension
    if numel(inter.lind_r1_rates)~=numel(sys.isotopes)
        error('the number of elements in inter.lind_r1_rates must be equal to the number of spins.');
    end
    
    % Enforce Lindblad theory if inter.lind_r1_rates rates are specified
    if ~ismember('lindblad',inter.relaxation)
        error('inter.lind_r1_rates can only be specified with Lindblad relaxation theory.');
    end
    
end

% Check R2 rates for T1/T2
if isfield(inter,'r2_rates')
    
    % Check type and values
    if (~isvector(inter.r2_rates))||(~isnumeric(inter.r2_rates))||...
         any(inter.r2_rates<0)||any(~isreal(inter.r2_rates))
        error('inter.r2_rates must be a vector of positive real numbers.');
    end
    
    % Check dimension
    if numel(inter.r2_rates)~=numel(sys.isotopes)
        error('the number of elements in inter.r2_rates must be equal to the number of spins.');
    end
    
    % Enforce T1,T2 theory if inter.r2_rates rates are specified
    if ~ismember('t1_t2',inter.relaxation)
        error('inter.r2_rates can only be specified with T1,T2 relaxation theory.');
    end
    
end

% Check R2 rates for Lindblad
if isfield(inter,'lind_r2_rates')
    
    % Check type and values
    if (~isvector(inter.lind_r2_rates))||(~isnumeric(inter.lind_r2_rates))||...
         any(inter.lind_r2_rates<0)||any(~isreal(inter.lind_r2_rates))
        error('inter.lind_r2_rates must be a vector of positive real numbers.');
    end
    
    % Check dimension
    if numel(inter.lind_r2_rates)~=numel(sys.isotopes)
        error('the number of elements in inter.lind_r2_rates must be equal to the number of spins.');
    end
    
    % Enforce Lindblad theory if inter.lind_r2_rates rates are specified
    if ~ismember('lindblad',inter.relaxation)
        error('inter.lind_r2_rates can only be specified with Lindblad relaxation theory.');
    end
    
end

% Check R1e rate for Weizmann DNP theory
if isfield(inter,'weiz_r1e')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r1e))||any(inter.weiz_r1e<0)||any(~isreal(inter.weiz_r1e))
        error('inter.weiz_r1e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r1e)~=1
        error('inter.weiz_r1e must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r1e is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r1e can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R1e rate for Nottingham DNP theory
if isfield(inter,'nott_r1e')
    
    % Check type and value
    if (~isnumeric(inter.nott_r1e))||any(inter.nott_r1e<0)||any(~isreal(inter.nott_r1e))
        error('inter.nott_r1e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r1e)~=1
        error('inter.nott_r1e must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1e is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r1e can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R2e rate for Weizmann DNP theory
if isfield(inter,'weiz_r2e')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r2e))||any(inter.weiz_r2e<0)||any(~isreal(inter.weiz_r2e))
        error('inter.weiz_r2e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r2e)~=1
        error('inter.weiz_r1e must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r2e is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r2e can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R2e rate for Nottingham DNP theory
if isfield(inter,'nott_r2e')
    
    % Check type and value
    if (~isnumeric(inter.nott_r2e))||any(inter.nott_r2e<0)||any(~isreal(inter.nott_r2e))
        error('inter.nott_r2e must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r2e)~=1
        error('inter.nott_r2e must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1e is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r2e can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R1n rate for Weizmann DNP theory
if isfield(inter,'weiz_r1n')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r1n))||any(inter.weiz_r1n<0)||any(~isreal(inter.weiz_r1n))
        error('inter.weiz_r1n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r1n)~=1
        error('inter.weiz_r1n must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r1n is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r1n can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R1n rate for Nottingham DNP theory
if isfield(inter,'nott_r1n')
    
    % Check type and value
    if (~isnumeric(inter.nott_r1n))||any(inter.nott_r1n<0)||any(~isreal(inter.nott_r1n))
        error('inter.nott_r1n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r1n)~=1
        error('inter.nott_r1n must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1n is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r1n can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R2n rate for Weizmann DNP theory
if isfield(inter,'weiz_r2n')
    
    % Check type and value
    if (~isnumeric(inter.weiz_r2n))||any(inter.weiz_r2n<0)||any(~isreal(inter.weiz_r2n))
        error('inter.weiz_r2n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.weiz_r2n)~=1
        error('inter.weiz_r1e must have a single element.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r2n is specified
    if ~ismember('weizmann',inter.relaxation)
        error('inter.weiz_r2n can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R2n rate for Nottingham DNP theory
if isfield(inter,'nott_r2n')
    
    % Check type and value
    if (~isnumeric(inter.nott_r2n))||any(inter.nott_r2n<0)||any(~isreal(inter.nott_r2n))
        error('inter.nott_r2n must be a positive real number.');
    end
    
    % Check dimension
    if numel(inter.nott_r2n)~=1
        error('inter.nott_r2n must have a single element.');
    end
    
    % Enforce Nottingham theory if inter.weiz_r1e is specified
    if ~ismember('nottingham',inter.relaxation)
        error('inter.nott_r2n can only be specified with Nottingham DNP relaxation theory.');
    end
    
end

% Check R1d rates for Weizmann DNP theory
if isfield(inter,'weiz_r1d')
    
    % Check type and values
    if (~isnumeric(inter.weiz_r1d))||any(inter.weiz_r1d(:)<0)||(~isreal(inter.weiz_r1d))
        error('inter.weiz_r1d must be a matrix of non-negative real numbers.');
    end
    
    % Check dimension
    if any(size(inter.weiz_r1d)~=[numel(sys.isotopes) numel(sys.isotopes)])
        error('both dimensions of inter.weiz_r1d must be equal to the number of spins.');
    end
    
    % Enforce Weizmann theory if inter.weiz_r1d rates are specified
    if ~strcmp(inter.relaxation,'weizmann')
        error('inter.weiz_r1d can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check R2d rates
if isfield(inter,'weiz_r2d')
    
    % Check type and values
    if (~isnumeric(inter.weiz_r2d))||any(inter.weiz_r2d(:)<0)||(~isreal(inter.weiz_r2d))
        error('inter.weiz_r2d must be a matrix of non-negative real numbers.');
    end
    
    % Check dimension
    if any(size(inter.weiz_r2d)~=[numel(sys.isotopes) numel(sys.isotopes)])
        error('both dimensions of inter.weiz_r2d must be equal to the number of spins.');
    end
    
    % Enforce Weizmann theory if R2d rates are specified
    if ~strcmp(inter.relaxation,'weizmann')
        error('inter.weiz_r2d can only be specified with Weizmann DNP relaxation theory.');
    end
    
end

% Check chemical kinetics
if isfield(inter,'chem')
    
    % Check reaction specifications
    if isfield(inter.chem,'parts')&&(~isfield(inter.chem,'rates'))
        error('reaction rates (inter.chem.rates) must be provided.');
    elseif isfield(inter.chem,'rates')&&(~isfield(inter.chem,'parts'))
        error('subsystem identifiers (inter.chem.parts) must be provided.');
    elseif isfield(inter.chem,'rates')&&(~isfield(inter.chem,'concs'))
        error('initial concentrations (inter.chem.concs) must be provided.');
    elseif isfield(inter.chem,'concs')&&(~isfield(inter.chem,'rates'))
        error('reaction rates (inter.chem.rates) must be provided.');
    end
    
    % Check chemical species specification
    if isfield(inter.chem,'parts')
        if ~iscell(inter.chem.parts)||(~all(cellfun(@isvector,inter.chem.parts)))
            error('inter.chem.parts must be a cell array of vectors.');
        end
        for n=1:numel(inter.chem.parts)
            for k=1:numel(inter.chem.parts)
                if (n~=k)&&(~isempty(intersect(inter.chem.parts{n},inter.chem.parts{k})))
                    error('a given spin can only belong to one chemical subsystem in inter.chem.parts variable.');
                end
            end
            if any(inter.chem.parts{n}<1)||any(inter.chem.parts{n}>numel(sys.isotopes))||...
               any(mod(inter.chem.parts{n},1)~=0)||(numel(unique(inter.chem.parts{n}))~=numel(inter.chem.parts{n}))
                error('elements of inter.chem.parts must be vectors of unique positive integers not exceeding the total number of spins.');
            end
            for k=1:numel(inter.chem.parts)
                if numel(inter.chem.parts{n})~=numel(inter.chem.parts{k})
                    error('all chemical subsystems must have the same number of spins.');
                end
            end
            for k=1:numel(inter.chem.parts)
                for m=1:numel(inter.chem.parts{n})
                    if ~strcmp(sys.isotopes{inter.chem.parts{n}(m)},sys.isotopes{inter.chem.parts{k}(m)})
                        error('isotope sequences in all chemical subsystems must be the same.');
                    end
                end
            end
        end
    end
    
    % Check reaction rate matrix
    if isfield(inter.chem,'rates')
        if (~isnumeric(inter.chem.rates))||(~isreal(inter.chem.rates))||(size(inter.chem.rates,1)~=size(inter.chem.rates,2))
            error('inter.chem.rates must be a real square matrix.');
        end
        if any(size(inter.chem.rates)~=numel(inter.chem.parts))
            error('both dimensions of inter.chem.rates matrix must be equal to the number of chemical subsystems.');
        end
    end
    
    % Check initial concentrations
    if isfield(inter.chem,'concs')
        if (~isnumeric(inter.chem.concs))||(~isreal(inter.chem.concs))||any(inter.chem.concs(:)<0)
            error('inter.chem.concs must be a vector of non-negative real numbers.');
        end
        if numel(inter.chem.concs)~=numel(inter.chem.parts)
            error('the number of initial concentrations must be equal to the number of chemical species.');
        end
    end
    
    % Check flux specifications
    if isfield(inter.chem,'flux_rate')&&(~isfield(inter.chem,'flux_type'))
        error('flux type (inter.chem.flux_type) must be provided.');
    elseif isfield(inter.chem,'flux_type')&&(~isfield(inter.chem,'flux_rate'))
        error('flux rates (inter.chem.flux_rate) must be provided.');
    end
    
    % Check flux rate matrix
    if isfield(inter.chem,'flux_rate')
        if (~isnumeric(inter.chem.flux_rate))||(~isreal(inter.chem.flux_rate))||...
           (size(inter.chem.flux_rate,1)~=size(inter.chem.flux_rate,2))
            error('inter.chem.flux_rate must be a real square matrix.');
        end
        if any(size(inter.chem.flux_rate)~=numel(sys.isotopes))
            error('both dimensions of inter.chem.flux_rate matrix must be equal to the number of spins.');
        end
    end
    
    % Check flux type
    if isfield(inter.chem,'flux_type')
        if ~ischar(inter.chem.flux_type)
            error('inter.chem.flux_type must be a character string.');
        end
        if ~ismember(inter.chem.flux_type,{'intermolecular','intramolecular'})
            error('incorrect flux type specification.');
        end
    end
    
    % Check radical pair kinetics
    if isfield(inter.chem,'rp_theory')
        if ~ischar(inter.chem.rp_theory)
            error('inter.chem.rp_theory must be a string.');
        end
        if ~ismember(inter.chem.rp_theory,{'haberkorn','jones-hore','exponential'})
            error('allowed values for inter.chem.rp_theory are ''exponential'', ''haberkorn'' and ''jones-hore''.');
        end
        if (~isfield(inter.chem,'rp_electrons'))||(~isfield(inter.chem,'rp_rates'))
            error('inter.chem.rp_electrons and inter.chem.rp_rates must be specified alongside inter.chem.rp_theory parameter.');
        end
    end
    if isfield(inter.chem,'rp_electrons')
        if (~isfield(inter.chem,'rp_theory'))||(~isfield(inter.chem,'rp_rates'))
            error('inter.chem.rp_theory and inter.chem.rp_rates must be specified alongside inter.chem.rp_electrons parameter.');
        end
        if (~isnumeric(inter.chem.rp_electrons))||(numel(inter.chem.rp_electrons)~=2)||any(ceil(inter.chem.rp_electrons)~=floor(inter.chem.rp_electrons))||any(inter.chem.rp_electrons<1)
            error('inter.chem.rp_electrons must be a vector of two positive integers.');
        end
        if any(inter.chem.rp_electrons>numel(sys.isotopes))||any(~cellfun(@(x)strcmp(x(1),'E'),sys.isotopes(inter.chem.rp_electrons)))
            error('at least one of the elements of inter.chem.rp_electrons does not refer to an electron.');
        end
    end
    if isfield(inter.chem,'rp_rates')
        if (~isfield(inter.chem,'rp_theory'))||(~isfield(inter.chem,'rp_electrons'))
            error('inter.chem.rp_theory and inter.chem.rp_electrons must be specified alongside inter.chem.rp_rates parameter.');
        end
        if (~isnumeric(inter.chem.rp_rates))||(numel(inter.chem.rp_rates)~=2)||(~isreal(inter.chem.rp_rates))||any(inter.chem.rp_rates<0)
            error('inter.chem.rp_rates must be a vector of two non-negative real numbers.');
        end
    end
    
end

end

% Those who beat their swords into plowshares will till 
% the soil for those who did not.
%
% Benjamin Franklin

