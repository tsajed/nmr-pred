% Tolerances and fundamental constants. Sets the various accuracy cut-offs,
% constants and tolerances used throughout Spinach kernel. 
%
% Modifications to this function are discouraged -- the accuracy settings
% should be modified by setting the sub-fields of the sys.tols structure,
% see the input preparation manual in the /docs directory.
%
% i.kuprov@soton.ac.uk
% hannah.hogben@chem.ox.ac.uk

function spin_system=tolerances(spin_system,sys)

% Interaction clean-up tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'inter_cutoff')
    spin_system.tols.inter_cutoff=sys.tols.inter_cutoff;
    report(spin_system,[pad('Drop interaction tensors with norms below (rad/s)',80) pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.inter_cutoff=eps();
    report(spin_system,[pad('Drop interaction tensors with norms below (rad/s)',80) pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (paranoid)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.inter_cutoff=1e-2;
    report(spin_system,[pad('Drop interaction tensors with norms below (rad/s)',80) pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.inter_cutoff=1e-10;
    report(spin_system,[pad('Drop interaction tensors with norms below (rad/s)',80) pad(num2str(spin_system.tols.inter_cutoff,'%0.8e'),20) ' (safe default)']);
end

% Liouvillian zero tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'liouv_zero')
    spin_system.tols.liouv_zero=sys.tols.liouv_zero;
    report(spin_system,[pad('Drop Liovillian terms with amplitudes below (rad/s)',80) pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.liouv_zero=eps();
    report(spin_system,[pad('Drop Liovillian terms with amplitudes below (rad/s)',80) pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (paranoid)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.liouv_zero=1e-5;
    report(spin_system,[pad('Drop Liovillian terms with amplitudes below (rad/s)',80) pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.liouv_zero=1e-10;
    report(spin_system,[pad('Drop Liovillian terms with amplitudes below (rad/s)',80) pad(num2str(spin_system.tols.liouv_zero,'%0.8e'),20) ' (safe default)']);
end

% Zero tolerance for the series terms in the exponential propagator
if isfield(sys,'tols')&&isfield(sys.tols,'prop_chop')
    spin_system.tols.prop_chop=sys.tols.prop_chop;
    report(spin_system,[pad('Drop exponential propagator terms with amplitudes below',80) pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.prop_chop=eps();
    report(spin_system,[pad('Drop exponential propagator terms with amplitudes below',80) pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (paranoid)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.prop_chop=1e-8;
    report(spin_system,[pad('Drop exponential propagator terms with amplitudes below',80) pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.prop_chop=1e-10;
    report(spin_system,[pad('Drop exponential propagator terms with amplitudes below',80) pad(num2str(spin_system.tols.prop_chop,'%0.8e'),20) ' (safe default)']);
end

% Subspace drop population tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'subs_drop')
    spin_system.tols.subs_drop=sys.tols.subs_drop;
    report(spin_system,[pad('Drop subspaces with total populations below',80) pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.subs_drop=eps();
    report(spin_system,[pad('Drop subspaces with total populations below',80) pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (paranoid)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.subs_drop=1e-2;
    report(spin_system,[pad('Drop subspaces with total populations below',80) pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.subs_drop=1e-10;
    report(spin_system,[pad('Drop subspaces with total populations below',80) pad(num2str(spin_system.tols.subs_drop,'%0.8e'),20) ' (safe default)']);
end

% Irrep population tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'irrep_drop')
    spin_system.tols.irrep_drop=sys.tols.irrep_drop;
    report(spin_system,[pad('Drop irreducible representations with populations below',80) pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.irrep_drop=eps();
    report(spin_system,[pad('Drop irreducible representations with populations below',80) pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (paranoid)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.irrep_drop=1e-2;
    report(spin_system,[pad('Drop irreducible representations with populations below',80) pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.irrep_drop=1e-10;
    report(spin_system,[pad('Drop irreducible representations with populations below',80) pad(num2str(spin_system.tols.irrep_drop,'%0.8e'),20) ' (safe default)']);
end

% Path drop tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'path_drop')
    spin_system.tols.path_drop=sys.tols.path_drop;
    report(spin_system,[pad('Disconect subspaces with cross-terms below (rad/s)',80) pad(num2str(spin_system.tols.path_drop,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.path_drop=eps();
    report(spin_system,[pad('Disconect subspaces with cross-terms below (rad/s)',80) pad(num2str(spin_system.tols.path_drop,'%0.8e'),20) ' (paranoid)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.path_drop=1e-2;
    report(spin_system,[pad('Disconect subspaces with cross-terms below (rad/s)',80) pad(num2str(spin_system.tols.path_drop,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.path_drop=1e-10;
    report(spin_system,[pad('Disconect subspaces with cross-terms below (rad/s)',80) pad(num2str(spin_system.tols.path_drop,'%0.8e'),20) ' (safe default)']);
end

% ZTE sample length
if isfield(sys,'tols')&&isfield(sys.tols,'zte_nsteps')
    spin_system.tols.zte_nsteps=sys.tols.zte_nsteps;
    report(spin_system,[pad('Number of steps in the zero track elimination sample',80) pad(num2str(spin_system.tols.zte_nsteps),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.zte_nsteps=16;
    report(spin_system,[pad('Number of steps in the zero track elimination sample',80) pad(num2str(spin_system.tols.zte_nsteps),20) ' (loose)']);
else
    spin_system.tols.zte_nsteps=32;
    report(spin_system,[pad('Number of steps in the zero track elimination sample',80) pad(num2str(spin_system.tols.zte_nsteps),20) ' (safe default)']);
end

% ZTE zero track tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'zte_tol')
    spin_system.tols.zte_tol=sys.tols.zte_tol;
    report(spin_system,[pad('Consider trajectory tracks to be zero below',80) pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('cowboy',sys.enable)
    spin_system.tols.zte_tol=1e-6;
    report(spin_system,[pad('Consider trajectory tracks to be zero below',80) pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (loose)']);
else
    spin_system.tols.zte_tol=1e-24;
    report(spin_system,[pad('Consider trajectory tracks to be zero below',80) pad(num2str(spin_system.tols.zte_tol,'%0.8e'),20) ' (safe default)']);
end

% ZTE state vector density threshold 
if isfield(sys,'tols')&&isfield(sys.tols,'zte_maxden')
    spin_system.tols.zte_maxden=sys.tols.zte_maxden;
    report(spin_system,[pad('ZTE off for state vector densities above',80) pad(num2str(spin_system.tols.zte_maxden),20) ' (user-specified)']);
else
    spin_system.tols.zte_maxden=0.5;
    report(spin_system,[pad('ZTE off for state vector densities above',80) pad(num2str(spin_system.tols.zte_maxden),20) ' (safe default)']);
end

% Proximity tolerance for dipolar couplings
if isfield(sys,'tols')&&isfield(sys.tols,'prox_cutoff')
    spin_system.tols.prox_cutoff=sys.tols.prox_cutoff;
    report(spin_system,[pad('Drop dipolar couplings over distances greater than (Angstrom)',80) pad(num2str(spin_system.tols.prox_cutoff),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.prox_cutoff=inf();
    report(spin_system,[pad('Drop dipolar couplings over distances greater than (Angstrom)',80) pad(num2str(spin_system.tols.prox_cutoff),20) ' (paranoid)']);
else
    spin_system.tols.prox_cutoff=100;
    report(spin_system,[pad('Drop dipolar couplings over distances greater than (Angstrom)',80) pad(num2str(spin_system.tols.prox_cutoff),20) ' (safe default)']);
end

% Krylov method switchover
if isfield(sys,'tols')&&isfield(sys.tols,'krylov_switchover')
    spin_system.tols.krylov_switchover=sys.tols.krylov_switchover;
    report(spin_system,[pad('Use Krylov propagation for nnz(L) above',80) pad(num2str(spin_system.tols.krylov_switchover),20) ' (user-specified)']);    
else
    spin_system.tols.krylov_switchover=1e5;
    report(spin_system,[pad('Use Krylov propagation for nnz(L) above',80) pad(num2str(spin_system.tols.krylov_switchover),20) ' (safe default)']);
end

% Basis printing hush tolerance
if isfield(sys,'tols')&&isfield(sys.tols,'basis_hush')
    spin_system.tols.basis_hush=sys.tols.basis_hush;
    report(spin_system,[pad('Do not print the basis for state space dimensions over',80) pad(num2str(spin_system.tols.basis_hush),20) ' (user-specified)']);
else
    spin_system.tols.basis_hush=256;
    report(spin_system,[pad('Do not print the basis for state space dimensions over',80) pad(num2str(spin_system.tols.basis_hush),20) ' (safe default)']);
end

% Subspace bundle size
if isfield(sys,'tols')&&isfield(sys.tols,'merge_dim')
    spin_system.tols.merge_dim=sys.tols.merge_dim;
    report(spin_system,[pad('Collect small subspaces into bundles of dimension',80) pad(num2str(spin_system.tols.merge_dim),20) ' (user-specified)']);
else
    spin_system.tols.merge_dim=1000;
    report(spin_system,[pad('Collect small subspaces into bundles of dimension',80) pad(num2str(spin_system.tols.merge_dim),20) ' (safe default)']);
end

% Sparse algebra tolerance on density
if isfield(sys,'tols')&&isfield(sys.tols,'dense_matrix')
    spin_system.tols.dense_matrix=sys.tols.dense_matrix;
    report(spin_system,[pad('Force sparse algebra for matrix density below',80) pad(num2str(spin_system.tols.dense_matrix),20) ' (user-specified)']);
else
    spin_system.tols.dense_matrix=0.10;
    report(spin_system,[pad('Force sparse algebra for matrix density below',80) pad(num2str(spin_system.tols.dense_matrix),20) ' (safe default)']);
end

% Sparse algebra tolerance on dimension
if isfield(sys,'tols')&&isfield(sys.tols,'small_matrix')
    spin_system.tols.small_matrix=sys.tols.small_matrix;
    report(spin_system,[pad('Force full algebra for matrix dimension below',80) pad(num2str(spin_system.tols.small_matrix),20) ' (user-specified)']);
else
    spin_system.tols.small_matrix=200;
    report(spin_system,[pad('Force full algebra for matrix dimension below',80) pad(num2str(spin_system.tols.small_matrix),20) ' (safe default)']);
end

% Relative accuracy of the elements of Redfield superoperator
if isfield(sys,'tols')&&isfield(sys.tols,'rlx_integration')
    spin_system.tols.rlx_integration=sys.tols.rlx_integration;
    report(spin_system,[pad('Relative accuracy for the elements of Redfield superoperator',80) pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (user-specified)']);
elseif isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.tols.rlx_integration=1e-6;
    report(spin_system,[pad('Relative accuracy for the elements of Redfield superoperator',80) pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (paranoid)']);
else
    spin_system.tols.rlx_integration=1e-4;
    report(spin_system,[pad('Relative accuracy for the elements of Redfield superoperator',80) pad(num2str(spin_system.tols.rlx_integration,'%0.8e'),20) ' (safe default)']);
end

% Algorithm selection for propagator derivatives
if isfield(sys,'tols')&&isfield(sys.tols,'dP_method')
    spin_system.tols.dP_method=sys.tols.dP_method;
    report(spin_system,[pad('Matrix exponential differentiation algorithm',80) pad(spin_system.tols.dP_method,20) ' (user-specified)']);
else
    spin_system.tols.dP_method='auxmat';
    report(spin_system,[pad('Matrix exponential differentiation algorithm',80) pad(spin_system.tols.dP_method,20) ' (safe default)']);
end

% Number of PBC images for dipolar couplings
if isfield(sys,'tols')&&isfield(sys.tols,'dd_ncells')
    spin_system.tols.dd_ncells=sys.tols.dd_ncells;
    report(spin_system,[pad('Number of PBC images for dipolar couplings in periodic systems',80) pad(num2str(spin_system.tols.dd_ncells),20) ' (user-specified)']);
else
    spin_system.tols.dd_ncells=2;
    report(spin_system,[pad('Number of PBC images for dipolar couplings in periodic systems',80) pad(num2str(spin_system.tols.dd_ncells),20) ' (safe default)']);
end

% Fundamental constants
spin_system.tols.hbar=1.054571628e-34;
report(spin_system,[pad('Planck constant (hbar)',80) pad(num2str(spin_system.tols.hbar,'%0.8e'),20)]);
spin_system.tols.kbol=1.3806503e-23;
report(spin_system,[pad('Boltzmann constant (k)',80) pad(num2str(spin_system.tols.kbol,'%0.8e'),20)]);
spin_system.tols.freeg=2.0023193043622;
report(spin_system,[pad('Free electron g-factor',80) pad(num2str(spin_system.tols.freeg,'%0.8e'),20)]);
spin_system.tols.mu0=4*pi*1e-7;
report(spin_system,[pad('Vacuum permeability',80) pad(num2str(spin_system.tols.mu0,'%0.8e'),20)]);

% Paranoia switches
if isfield(sys,'enable')&&ismember('paranoia',sys.enable)
    spin_system.sys.disable=unique([spin_system.sys.disable {'trajlevel','krylov','clean-up','expv'}]);
end

end

% Man once surrendering his reason, has no remaining guard against
% absurdities the most monstrous, and like a ship without rudder, 
% is the sport of every wind. 
%
% Thomas Jefferson

