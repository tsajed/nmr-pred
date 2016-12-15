% Kernel integrity control. Checks for collisions between Spinach func-
% tions and anything else that the user may have installed or written
% in the current Matlab instance.
%
% Do not switch this off -- collisions of function names and path probl-
% ems are the most frequent support topic at the forum.
%
% i.kuprov@soton.ac.uk
% dg1g13@soton.ac.uk (David Goodwin)

function existentials()

% Do not run existential testing inside parallel pools
if isworkernode, return; end

% Locate Spinach
root_dir=which('existentials'); root_dir=root_dir(1:end-32);

% Check location of kernel functions
kernel_functions={'assume','average','basis','coherence','correlation','create',...
                  'carrier','decouple','equilibrium','evolution','execute',...
                  'hamiltonian','homospoil','kinetics','frqoffset','krylov',...
                  'lindbladian','operator','orientation','propagator','rotframe',...
                  'reduce','relaxation','residual','spin','splice','state',...
                  'step','trajan','trajsimil','spinlock','stateinfo','stepsize',...
                  'thermalize','singlet'};
for n=1:numel(kernel_functions)
    expected_location=[root_dir filesep 'kernel' filesep kernel_functions{n} '.m'];
    apparent_location=which(kernel_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' kernel_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of kernel utilities
kernel_utilities={'absorb','anax2dcm','anax2quat','ang2cgsppm','apodization','axis_1d','axrh2mat','banner',...
                  'binpack','castep2nqi','cce','centroid','cgsppm2ang','clean_up','clebsch_gordan',...
                  'conmat','corrfun','crop_2d','dcm2euler','dcm2quat','dcm2wigner','dfpt','dihedral',...
                  'dilute','dipolar','dirdiff','eeqq2nqi','equilibrate','euler2dcm','euler2wigner',...
                  'existentials','expmint','fdhess','fdkup','fdlap','fdmat','fdvec','fdweights','fftdiff',...
                  'fourlap','fpl2phan','fpl2rho','fprint','frac2cart','gauss2mhz','gaussfun','hartree2joule',...
                  'hdot','hess_reorder','hilb2liouv','human2opspec','int_2d','interpmat','intrep',...
                  'irr_sph_ten','isworkernode','kill_spin','killcross','killdiag','krondelta','lcurve',...
                  'lin2lm','lin2lmn','lm2lin','lmn2lin','lorentzfun','mat2sphten','md5_hash','mhz2gauss',...
                  'molplot','mprealloc','mt2hz','negligible','p_superop','pad','path_trace','pauli',...
                  'perm_group','phan2fpl','plot_1d','plot_2d','plot_3d','probmax','quat2anax','quat2dcm',...
                  'relax_split','report','rotor_stack','scomponents','shift_iso','significant','slice_2d',...
                  'sniff','sparse2csr','spher_harmon','sphten2mat','sphten2oper','sphten2state','sphten2zeeman',...
                  'statmerge','stitch','summary','svd_shrink','sweep2ticks','symmetrize','symmetry',...
                  'tensor_analysis','tolerances','transfermat','unit_oper','unit_state','v2fplanck',...
                  'volplot','wigner','xyz2dd','xyz2hfc','xyz2sph','zeeman2sphten','zfs2mat','zte'};
for n=1:numel(kernel_utilities)
    expected_location=[root_dir filesep 'kernel' filesep 'utilities' filesep kernel_utilities{n} '.m'];
    apparent_location=which(kernel_utilities{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' kernel_utilities{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of pulse shapes
kernel_shapes={'cartesian2polar','chirp_pulse_af','chirp_pulse_xy','grad_pulse','grad_sandw','polar2cartesian',...
               'pulse_shape','read_wave','sawtooth','shaped_pulse_af','shaped_pulse_xy','sim_pulse','triwave',...
               'vg_pulse','wave_basis'};
for n=1:numel(kernel_shapes)
    expected_location=[root_dir filesep 'kernel' filesep 'pulses' filesep kernel_shapes{n} '.m'];
    apparent_location=which(kernel_shapes{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' kernel_shapes{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of optimal control functions
oc_functions={'control_sys','fminkrotov','fminnewton','fminsimplex','grape','hessprep','hessreg','krotov',...
              'lbfgs','linesearch','optim_report','optim_tols','penalty','quasinewton'};
for n=1:numel(oc_functions)
    expected_location=[root_dir filesep 'kernel' filesep 'optimcon' filesep oc_functions{n} '.m'];
    apparent_location=which(oc_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' oc_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of grid functions
grid_functions={'gaussleg','grid_kron','grid_test','repulsion','shrewd'};
for n=1:numel(grid_functions)
    expected_location=[root_dir filesep 'kernel' filesep 'grids' filesep grid_functions{n} '.m'];
    apparent_location=which(grid_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' grid_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of external functions
ext_functions={'assofix','b2r','combnk','expv','fourdif','jacobianest','lgwt','lpredict','phantom3d','simps'};
for n=1:numel(ext_functions)
    expected_location=[root_dir filesep 'kernel' filesep 'external' filesep ext_functions{n} '.m'];
    apparent_location=which(ext_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' ext_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of cache functions
cache_functions={'ist_product_table','sle_operators'};
for n=1:numel(cache_functions)
    expected_location=[root_dir filesep 'kernel' filesep 'cache' filesep cache_functions{n} '.m'];
    apparent_location=which(cache_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' cache_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of context functions
context_functions={'crystal','doublerot','floquet','gridfree','hydrodynamics','imaging','liquid',...
                   'oscillator','powder','roadmap','singlerot'};
for n=1:numel(cache_functions)
    expected_location=[root_dir filesep 'kernel' filesep 'contexts' filesep context_functions{n} '.m'];
    apparent_location=which(context_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' context_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of experiments functions
exp_functions={['hyperpol' filesep 'dnp_field_scan'],['hyperpol' filesep 'dnp_freq_scan'],...
               ['hyperpol' filesep 'masdnp'],['hyperpol' filesep 'solid_effect'],...
               ['overtone' filesep 'overtone_a'],['overtone' filesep 'overtone_cp'],...
               ['overtone' filesep 'overtone_dante'],['overtone' filesep 'overtone_pa'],...
               ['pseudocon' filesep 'cg_fast'],['pseudocon' filesep 'chi_eff'],...
               ['pseudocon' filesep 'hfc2pcs'],['pseudocon' filesep 'hfc2pms'],...
               ['pseudocon' filesep 'ilpcs'],['pseudocon' filesep 'ipcs'],...
               ['pseudocon' filesep 'ippcs'],['pseudocon' filesep 'kpcs'],...
               ['pseudocon' filesep 'logfactorial'],['pseudocon' filesep 'lpcs'],...
               ['pseudocon' filesep 'pcs_combi_fit'],['pseudocon' filesep 'pcs2chi'],...
               ['pseudocon' filesep 'pms2chi'],['pseudocon' filesep 'points2mult'],...
               ['pseudocon' filesep 'ppcs'],['pseudocon' filesep 'qform2sph'],...
               'acquire','clip_hsqc','cn2d',...
               'cosy','crazed','crosspol','dante','deer_3p_hard_deer','deer_3p_hard_echo','deer_3p_soft_deer',...
               'deer_3p_soft_diag','deer_3p_soft_hole','deer_4p_soft_deer','deer_4p_soft_diag','deer_4p_soft_hole',...
               'deer_analyt','dqf_cosy','endor_cw','endor_davies','endor_mims','eseem','fieldsweep','hetcor',...
               'hmqc','hnco','hncoca','holeburn','hp_acquire','hsqc','hyscore','lcosy','levelpop','m2s',...
               'noesy','noesyhsqc','oopeseem','rapidscan','roesy','rydmr','rydmr_exp','s2m','slowpass',...
               'sp_acquire','traject','ufcosy','zerofield'};
for n=1:numel(exp_functions)
    expected_location=[root_dir filesep 'experiments' filesep exp_functions{n} '.m'];
    apparent_location=which(exp_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' exp_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

% Check location of miscellaneous functions
etc_functions={'cst_display','destreak','fid2ascii','g2spinach','gparse',...
               'guess_csa_pro','guess_j_pro','guess_j_nuc','hfc_display','karplus_fit',...
               'nuclacid','protein','read_bmrb','read_pdb_pro','read_pdb_nuc','s2spinach',...
               'strychnine','zfs_sampling'};
for n=1:numel(etc_functions)
    expected_location=[root_dir filesep 'etc' filesep etc_functions{n} '.m'];
    apparent_location=which(etc_functions{n});
    if ~strcmp(expected_location,apparent_location)
        disp(['Function ' etc_functions{n}]);
        disp(['Expected location ' expected_location]);
        disp(['Apparent location ' apparent_location]);
        error(['function name collision: either Matlab path is not set correctly or you '...
               'have a file with the same name as one of the essential Spinach functions '...
               'somewhere, and it prevents Spinach from running.']);
    end
end

end

% It has been my observation that most people get ahead
% during the time that others waste.
%
% Henry Ford

