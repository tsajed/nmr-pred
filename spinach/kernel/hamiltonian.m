% Hamiltonian operator / superoperator and its rotational decomposition.
% Descriptor and operator generation is parallelized. Syntax:
%
%                [H,Q]=hamiltonian(spin_system,operator_type)
%
% In Liouville space calculations, operator_type can be set to:
%                      
%             'left' - produces left side product superoperator
%
%            'right' - produces right side product superoperator
%
%             'comm' - produces commutation superoperator (default)
%
%            'acomm' - produces anticommutation superoperator
%
% In Hilbert space calculations operator_type parameter is ignored.
%
% Outputs:
%                     H   - isotropic part
%
%                     Q   - twenty-five matrices giving irreducible
%                           components of the anisotropic part
%
% Note: this function returns the full rotational basis that contains 
%       information about all orientations. Use orientation.m to gene-
%       rate the Hamiltonian at a specific orientation. 
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk
% hannah.hogben@chem.ox.ac.uk
% d.savostyanov@soton.ac.uk

function [H,Q]=hamiltonian(spin_system,operator_type)

% Check consistency
grumble(spin_system);

% Set the default for the type
if ~exist('operator_type','var'), operator_type='comm'; end

% Decide if the Q part is required
build_aniso=(nargout>1);

% Inform the user
report(spin_system,'building Hamiltonian descriptor...');

% Preallocate spin indices
spin_L=zeros(spin_system.comp.nspins,8);
spin_S=zeros(spin_system.comp.nspins,8);

% Preallocate operator specifications
oper_L(1:spin_system.comp.nspins,1:8)={'E'};
oper_S(1:spin_system.comp.nspins,1:8)={'E'};

% Preallocate isotropic Hamiltonian coefficients
H_coeff=zeros(spin_system.comp.nspins,8);

% Preallocate spherical tensor coefficients
T_coeff=zeros(spin_system.comp.nspins,8,5);

% Preallocate ireducible components
irrcomp=zeros(spin_system.comp.nspins,8,5);

% Process Zeeman interactions
for n=1:spin_system.comp.nspins
    
    % Write the isotropic part
    switch spin_system.inter.zeeman.strength{n}
        
        case {'full','z_full'}
            
            % Add the carrier frequency
            zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3+...
                             spin_system.inter.basefrqs(n);
            
            % Update the Hamiltonian
            if significant(zeeman_iso,spin_system.tols.inter_cutoff)
                
                % Inform the user
                report(spin_system,['complete isotropic Zeeman interaction for spin ' num2str(n) '...']);
                report(spin_system,['           (Lz) x ' num2str(zeeman_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                spin_L(n,2)=n; oper_L(n,2)={'Lz'}; H_coeff(n,2)=zeeman_iso;
                
            end
            
        case {'secular','z_offs'}
            
            % Skip the carrier frequency
            zeeman_iso=trace(spin_system.inter.zeeman.matrix{n})/3;
            
            % Update the Hamiltonian
            if significant(zeeman_iso,spin_system.tols.inter_cutoff)
                
                % Inform the user
                report(spin_system,['offset isotropic Zeeman interaction for spin ' num2str(n) '...']);
                report(spin_system,['           (Lz) x ' num2str(zeeman_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                spin_L(n,2)=n; oper_L(n,2)={'Lz'}; H_coeff(n,2)=zeeman_iso;
                
            end
            
        case {'ignore','+','-'}
            
            % Inform the user
            report(spin_system,['isotropic Zeeman interaction ignored for spin ' num2str(n) '.']);
            
        otherwise
            
            % Bomb out with unexpected strength parameters
            error(['unknown strength specification for the Zeeman interaction of spin ' num2str(spin_L)]);
            
    end
    
    % Process anisotropic part if required
    if build_aniso
        
        % Get second rank spherical tensor components
        [~,~,phi_zeeman]=mat2sphten(spin_system.inter.zeeman.matrix{n});
        
        % Process Zeeman interactions
        if significant(phi_zeeman,spin_system.tols.inter_cutoff)
            
            % Store coefficients
            for k=1:3, irrcomp(n,k,:)=phi_zeeman; end
            
            % Write irreducible spherical tensors
            switch spin_system.inter.zeeman.strength{n}
                
                case 'full'
                    
                    % Inform the user
                    report(spin_system,['complete anisotropic Zeeman interaction for spin ' num2str(n) '...']);
                    report(spin_system,['                -0.5*(Lp) x ' num2str(phi_zeeman(2)/(2*pi))   ' Hz']);
                    report(spin_system,['           sqrt(2/3)*(Lz) x ' num2str(phi_zeeman(3)/(2*pi))   ' Hz']);
                    report(spin_system,['                 0.5*(Lm) x ' num2str(phi_zeeman(4)/(2*pi))   ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,1)=n; oper_L(n,1)={'L+'}; T_coeff(n,1,2)=-0.5;
                    spin_L(n,2)=n; oper_L(n,2)={'Lz'}; T_coeff(n,2,3)=sqrt(2/3);
                    spin_L(n,3)=n; oper_L(n,3)={'L-'}; T_coeff(n,3,4)=+0.5;
                    
                case {'secular','z_full','z_offs'}
                    
                    % Inform the user
                    report(spin_system,['Z part of the anisotropic Zeeman interaction for spin ' num2str(n) '...']);
                    report(spin_system,['           sqrt(2/3)*(Lz) x ' num2str(phi_zeeman(3)/(2*pi)) ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2)=n; oper_L(n,2)={'Lz'}; T_coeff(n,2,3)=sqrt(2/3);
                    
                    
                case '+'
                    
                    % Inform the user
                    report(spin_system,['"+" part of the anisotropic Zeeman interaction for spin ' num2str(n) '...']);
                    report(spin_system,['          -0.5*(Lp) x ' num2str(phi_zeeman(2)/(2*pi)) ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,1)=n; oper_L(n,1)={'L+'}; T_coeff(n,1,2)=-0.5;
                    
                case '-'
                    
                    % Inform the user
                    report(spin_system,['"-" part of the anisotropic Zeeman interaction for spin ' num2str(n) '...']);
                    report(spin_system,['           0.5*(Lm) x ' num2str(phi_zeeman(4)/(2*pi)) ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,3)=n; oper_L(n,3)={'L-'}; T_coeff(n,3,4)=+0.5;
                    
                case 'ignore'
                    
                    % Inform the user
                    report(spin_system,['anisotropic Zeeman interaction ignored for spin ' num2str(n) '.']);
                    
                otherwise
                    
                    % Bomb out with unexpected strength parameters
                    error(['unknown Zeeman interaction strength specification for spin ' num2str(n) '.']);
                    
            end
        end
        
        % Get second rank spherical tensor components
        [~,~,phi_quad]=mat2sphten(spin_system.inter.coupling.matrix{n,n});
        
        % Process quadrupolar interactions
        if significant(phi_quad,spin_system.tols.inter_cutoff)
            
            % Store coefficients
            for k=4:8, irrcomp(n,k,:)=phi_quad; end
            
            % Process the coupling
            switch spin_system.inter.coupling.strength{n,n}
                
                case 'strong'
                    
                    % Inform the user
                    report(spin_system,['complete quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2,+2) x ' num2str(phi_quad(1)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2,+1) x ' num2str(phi_quad(2)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2, 0) x ' num2str(phi_quad(3)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2,-1) x ' num2str(phi_quad(4)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           T(2,-2) x ' num2str(phi_quad(5)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
    
                    spin_L(n,4)=n; oper_L(n,4)={'T2,+2'}; T_coeff(n,4,1)=1;
                    spin_L(n,5)=n; oper_L(n,5)={'T2,+1'}; T_coeff(n,5,2)=1;
                    spin_L(n,6)=n; oper_L(n,6)={'T2,0'};  T_coeff(n,6,3)=1;
                    spin_L(n,7)=n; oper_L(n,7)={'T2,-1'}; T_coeff(n,7,4)=1;
                    spin_L(n,8)=n; oper_L(n,8)={'T2,-2'}; T_coeff(n,8,5)=1;
                    
                case {'secular','T2,0'}
                    
                    % Inform the user
                    report(spin_system,['secular part of the quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2, 0) x ' num2str(phi_quad(3)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,6)=n; oper_L(n,6)={'T2,0'};  T_coeff(n,6,3)=1;
                    
                case {'T2,+2'}
                    
                    % Inform the user
                    report(spin_system,['T(2,+2) part of the quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2,+2) x ' num2str(phi_quad(1)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,4)=n; oper_L(n,4)={'T2,+2'}; T_coeff(n,4,1)=1;
                    
                case {'T2,-2'}
                    
                    % Inform the user
                    report(spin_system,['T(2,-2) part of the quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2,-2) x ' num2str(phi_quad(5)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,8)=n; oper_L(n,8)={'T2,-2'}; T_coeff(n,8,5)=1;
                    
                case {'T2,+1'}
                    
                    % Inform the user
                    report(spin_system,['T(2,+1) part of the quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2,+1) x ' num2str(phi_quad(2)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,5)=n; oper_L(n,5)={'T2,+1'}; T_coeff(n,5,2)=1;
                    
                case {'T2,-1'}
                    
                    % Inform the user
                    report(spin_system,['T(2,-1) part of the quadratic coupling for spin ' num2str(n) '...']);
                    report(spin_system,['           T(2,-1) x ' num2str(phi_quad(4)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,7)=n; oper_L(n,7)={'T2,-1'}; T_coeff(n,7,4)=1;
                    
                case 'ignore'
                    
                    % Inform the user
                    report(spin_system,['quadratic coupling ignored for spin ' num2str(n) '.']);
                    
                otherwise
                    
                    % Bomb out with unexpected strength parameters
                    error(['unknown strength specification for the quadratic coupling of spin ' num2str(n)]);
                    
            end
            
        end
        
    end
    
end

% Pack single-spin descriptor table
D1=table(reshape(spin_L,[8*spin_system.comp.nspins 1]),reshape(spin_S,[8*spin_system.comp.nspins 1]),...
         reshape(oper_L,[8*spin_system.comp.nspins 1]),reshape(oper_S,[8*spin_system.comp.nspins 1]),...
         reshape(H_coeff,[8*spin_system.comp.nspins 1]),reshape(T_coeff,[8*spin_system.comp.nspins 5]),...
         reshape(irrcomp,[8*spin_system.comp.nspins 5]),'VariableNames',{'L','S','opL','opS','H','T','phi'});

% Kill insignificant rows
tol1=spin_system.tols.inter_cutoff; tol2=eps;
D1((abs(D1.H)<tol1)&((sum(abs(D1.phi),2)<tol1)|(sum(abs(D1.T),2)<tol2)),:)=[];

% Clean up variables for re-use
clear('spin_L','spin_S','oper_L','oper_S','H_coeff','T_coeff','irrcomp');

% Discover significant bilinear couplings and build interacting spin pair list
[L,S]=find(cellfun(@norm,spin_system.inter.coupling.matrix)>spin_system.tols.inter_cutoff);
quad_couplings=(L==S); L(quad_couplings)=[]; S(quad_couplings)=[];
pair_list=[L S]; if isempty(pair_list), pair_list=[]; end

% Preallocate spin indices
spin_L=zeros(size(pair_list,1),3,3);
spin_S=zeros(size(pair_list,1),3,3);

% Preallocate operator specifications
oper_L(1:size(pair_list,1),1:3,1:3)={'E'};
oper_S(1:size(pair_list,1),1:3,1:3)={'E'};

% Preallocate isotropic Hamiltonian coefficients
H_coeff=zeros(size(pair_list,1),3,3);

% Preallocate spherical tensor coefficients
T_coeff=zeros(size(pair_list,1),3,3,5);

% Preallocate ireducible components
irrcomp=zeros(size(pair_list,1),3,3,5);

% Loop over the pair list
for n=1:size(pair_list,1)
    
    % Extract spin numbers
    L=pair_list(n,1); S=pair_list(n,2);
    
    % Get the isotropic coupling constant
    coupling_iso=trace(spin_system.inter.coupling.matrix{L,S})/3;
    
    % Check if the coupling is significant
    if significant(coupling_iso,spin_system.tols.inter_cutoff)
        
        % Process the coupling
        switch spin_system.inter.coupling.strength{L,S}
            
            case {'strong','secular'}
                
                % Inform the user
                report(spin_system,['complete isotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                report(spin_system,['           (LxSx+LySy+LzSz) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; H_coeff(n,2,2)=coupling_iso;
                spin_L(n,1,3)=L; spin_S(n,1,3)=S; oper_L(n,1,3)={'L+'}; oper_S(n,1,3)={'L-'}; H_coeff(n,1,3)=coupling_iso/2;
                spin_L(n,3,1)=L; spin_S(n,3,1)=S; oper_L(n,3,1)={'L-'}; oper_S(n,3,1)={'L+'}; H_coeff(n,3,1)=coupling_iso/2;
                
            case {'weak','z*','*z','zz'}
                
                % Inform the user
                report(spin_system,['(z,z) part of the isotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                report(spin_system,['           (LzSz) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; H_coeff(n,2,2)=coupling_iso;
                
            case {'+-'}
                
                % Inform the user
                report(spin_system,['(+,-) part of the isotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                report(spin_system,['           0.5*(LpSm) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                spin_L(n,1,3)=L; spin_S(n,1,3)=S; oper_L(n,1,3)={'L+'}; oper_S(n,1,3)={'L-'}; H_coeff(n,1,3)=coupling_iso/2;
                
            case {'-+'}
                
                % Inform the user
                report(spin_system,['(-,+) part of the isotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                report(spin_system,['           0.5*(LmSp) x ' num2str(coupling_iso/(2*pi)) ' Hz']);
                
                % Update the Hamiltonian descriptor
                spin_L(n,3,1)=L; spin_S(n,3,1)=S; oper_L(n,3,1)={'L-'}; oper_S(n,3,1)={'L+'}; H_coeff(n,3,1)=coupling_iso/2;
                
            case {'ignore','z+','z-','+z','-z','++','--','T2,0','T(2,+1)','T(2,-1)','T(2,+2)','T(2,-2)'}
                
                % Inform the user
                report(spin_system,['isotropic coupling ignored for spins ' num2str(L) ',' num2str(S) '.']);
                
            otherwise
                
                % Bomb out with unexpected strength parameters
                error(['unknown strength specification for the bilinear coupling between spins ' num2str(L) ' and ' num2str(S)]);
        end
        
    end
    
    % Process anisotropic part if required
    if build_aniso
        
        % Get second rank spherical tensor components
        [~,~,phi_coupling]=mat2sphten(spin_system.inter.coupling.matrix{L,S});
        
        % Check if it is significant
        if significant(phi_coupling,spin_system.tols.inter_cutoff)
            
            % Store coefficients
            for k=1:3, for m=1:3, irrcomp(n,k,m,:)=phi_coupling; end; end
            
            % Write irreducible spherical tensors
            switch spin_system.inter.coupling.strength{L,S}
                
                case 'strong'
                    
                    % Inform the user
                    report(spin_system,['complete anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['                           0.5*(LpSp) x ' num2str(phi_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['                     -0.5*(LzSp+LpSz) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['    sqrt(2/3)*(LzSz-0.25*(LpSm+LmSp)) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['                      0.5*(LzSm+LmSz) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['                           0.5*(LmSm) x ' num2str(phi_coupling(5)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; T_coeff(n,2,2,3)=+sqrt(2/3);
                    spin_L(n,1,3)=L; spin_S(n,1,3)=S; oper_L(n,1,3)={'L+'}; oper_S(n,1,3)={'L-'}; T_coeff(n,1,3,3)=-sqrt(2/3)/4;
                    spin_L(n,3,1)=L; spin_S(n,3,1)=S; oper_L(n,3,1)={'L-'}; oper_S(n,3,1)={'L+'}; T_coeff(n,3,1,3)=-sqrt(2/3)/4;
                    spin_L(n,2,1)=L; spin_S(n,2,1)=S; oper_L(n,2,1)={'Lz'}; oper_S(n,2,1)={'L+'}; T_coeff(n,2,1,2)=-1/2;
                    spin_L(n,1,2)=L; spin_S(n,1,2)=S; oper_L(n,1,2)={'L+'}; oper_S(n,1,2)={'Lz'}; T_coeff(n,1,2,2)=-1/2;
                    spin_L(n,2,3)=L; spin_S(n,2,3)=S; oper_L(n,2,3)={'Lz'}; oper_S(n,2,3)={'L-'}; T_coeff(n,2,3,4)=+1/2;
                    spin_L(n,3,2)=L; spin_S(n,3,2)=S; oper_L(n,3,2)={'L-'}; oper_S(n,3,2)={'Lz'}; T_coeff(n,3,2,4)=+1/2;
                    spin_L(n,1,1)=L; spin_S(n,1,1)=S; oper_L(n,1,1)={'L+'}; oper_S(n,1,1)={'L+'}; T_coeff(n,1,1,1)=+1/2;
                    spin_L(n,3,3)=L; spin_S(n,3,3)=S; oper_L(n,3,3)={'L-'}; oper_S(n,3,3)={'L-'}; T_coeff(n,3,3,5)=+1/2;
                    
                case 'z*'
                    
                    % Inform the user
                    report(spin_system,['(z,*) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['                -0.5*(LzSp) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           sqrt(2/3)*(LzSz) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['                 0.5*(LzSm) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; T_coeff(n,2,2,3)=+sqrt(2/3);
                    spin_L(n,2,1)=L; spin_S(n,2,1)=S; oper_L(n,2,1)={'Lz'}; oper_S(n,2,1)={'L+'}; T_coeff(n,2,1,2)=-1/2;
                    spin_L(n,2,3)=L; spin_S(n,2,3)=S; oper_L(n,2,3)={'Lz'}; oper_S(n,2,3)={'L-'}; T_coeff(n,2,3,4)=+1/2;
                    
                case '*z'
                    
                    % Inform the user
                    report(spin_system,['(*,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['                -0.5*(LpSz) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['           sqrt(2/3)*(LzSz) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    report(spin_system,['                 0.5*(LmSz) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; T_coeff(n,2,2,3)=+sqrt(2/3);
                    spin_L(n,1,2)=L; spin_S(n,1,2)=S; oper_L(n,1,2)={'L+'}; oper_S(n,1,2)={'Lz'}; T_coeff(n,1,2,2)=-1/2;
                    spin_L(n,3,2)=L; spin_S(n,3,2)=S; oper_L(n,3,2)={'L-'}; oper_S(n,3,2)={'Lz'}; T_coeff(n,3,2,4)=+1/2;
                    
                case {'secular','T2,0'}
                    
                    % Inform the user
                    report(spin_system,['secular part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['   sqrt(2/3)*(LzSz-0.25*(LpSm+LmSp)) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; T_coeff(n,2,2,3)=+sqrt(2/3);
                    spin_L(n,1,3)=L; spin_S(n,1,3)=S; oper_L(n,1,3)={'L+'}; oper_S(n,1,3)={'L-'}; T_coeff(n,1,3,3)=-sqrt(2/3)/4;
                    spin_L(n,3,1)=L; spin_S(n,3,1)=S; oper_L(n,3,1)={'L-'}; oper_S(n,3,1)={'L+'}; T_coeff(n,3,1,3)=-sqrt(2/3)/4;
                    
                case {'weak','zz'}
                    
                    % Inform the user
                    report(spin_system,['(z,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['           sqrt(2/3)*(LzSz) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,2)=L; spin_S(n,2,2)=S; oper_L(n,2,2)={'Lz'}; oper_S(n,2,2)={'Lz'}; T_coeff(n,2,2,3)=+sqrt(2/3);
                    
                case 'z+'
                    
                    % Inform the user
                    report(spin_system,['(z,+) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['          -0.5*(LzSp) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,1)=L; spin_S(n,2,1)=S; oper_L(n,2,1)={'Lz'}; oper_S(n,2,1)={'L+'}; T_coeff(n,2,1,2)=-1/2;
                    
                case '+z'
                    
                    % Inform the user
                    report(spin_system,['(+,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['          -0.5*(LpSz) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,1,2)=L; spin_S(n,1,2)=S; oper_L(n,1,2)={'L+'}; oper_S(n,1,2)={'Lz'}; T_coeff(n,1,2,2)=-1/2;
                    
                case 'T(2,+1)'
                    
                    % Inform the user
                    report(spin_system,['T(2,+1) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['     -0.5*(LpSz+LzSp) x ' num2str(phi_coupling(2)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,1)=L; spin_S(n,2,1)=S; oper_L(n,2,1)={'Lz'}; oper_S(n,2,1)={'L+'}; T_coeff(n,2,1,2)=-1/2;
                    spin_L(n,1,2)=L; spin_S(n,1,2)=S; oper_L(n,1,2)={'L+'}; oper_S(n,1,2)={'Lz'}; T_coeff(n,1,2,2)=-1/2;
                    
                case 'z-'
                    
                    % Inform the user
                    report(spin_system,['(z,-) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['           0.5*(LzSm) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,3)=L; spin_S(n,2,3)=S; oper_L(n,2,3)={'Lz'}; oper_S(n,2,3)={'L-'}; T_coeff(n,2,3,4)=+1/2;
                    
                case '-z'
                    
                    % Inform the user
                    report(spin_system,['(-,z) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['           0.5*(LmSz) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,3,2)=L; spin_S(n,3,2)=S; oper_L(n,3,2)={'L-'}; oper_S(n,3,2)={'Lz'}; T_coeff(n,3,2,4)=+1/2;
                    
                case 'T(2,-1)'
                    
                    % Inform the user
                    report(spin_system,['T(2,-1) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['      0.5*(LmSz+LzSm) x ' num2str(phi_coupling(4)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,2,3)=L; spin_S(n,2,3)=S; oper_L(n,2,3)={'Lz'}; oper_S(n,2,3)={'L-'}; T_coeff(n,2,3,4)=+1/2;
                    spin_L(n,3,2)=L; spin_S(n,3,2)=S; oper_L(n,3,2)={'L-'}; oper_S(n,3,2)={'Lz'}; T_coeff(n,3,2,4)=+1/2;
                    
                case '+-'
                    
                    % Inform the user
                    report(spin_system,['(+,-) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['          -0.25*sqrt(2/3)*(LpSm) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,1,3)=L; spin_S(n,1,3)=S; oper_L(n,1,3)={'L+'}; oper_S(n,1,3)={'L-'}; T_coeff(n,1,3,3)=-sqrt(2/3)/4;
                    
                case '-+'
                    
                    % Inform the user
                    report(spin_system,['(-,+) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['          -0.25*sqrt(2/3)*(LmSp) x ' num2str(phi_coupling(3)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,3,1)=L; spin_S(n,3,1)=S; oper_L(n,3,1)={'L-'}; oper_S(n,3,1)={'L+'}; T_coeff(n,3,1,3)=-sqrt(2/3)/4;
                    
                case {'++','T2,+2'}
                    
                    % Inform the user
                    report(spin_system,['(+,+) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['           0.5*(LpSp) x ' num2str(phi_coupling(1)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,1,1)=L; spin_S(n,1,1)=S; oper_L(n,1,1)={'L+'}; oper_S(n,1,1)={'L+'}; T_coeff(n,1,1,1)=+1/2;
                    
                case {'--','T2,-2'}
                    
                    % Inform the user
                    report(spin_system,['(-,-) part of the anisotropic coupling for spins ' num2str(L) ',' num2str(S) '...']);
                    report(spin_system,['           0.5*(LmSm) x ' num2str(phi_coupling(5)/(2*pi),'%+0.5e') ' Hz']);
                    
                    % Prepare spherical tensor descriptors
                    spin_L(n,3,3)=L; spin_S(n,3,3)=S; oper_L(n,3,3)={'L-'}; oper_S(n,3,3)={'L-'}; T_coeff(n,3,3,5)=+1/2;
                    
                case 'ignore'
                    
                    % Inform the user
                    report(spin_system,['anisotropic coupling ignored for spins ' num2str(L) ',' num2str(S) '.']);
                    
                otherwise
                    
                    % Bomb out with unexpected strength parameters
                    error(['unknown strength specification for the bilinear coupling between spins ' num2str(L) ' and ' num2str(S)]);
                    
            end
        end
    end
end
    
% Pack two-spin descriptor table
D2=table(reshape(spin_L, [9*size(pair_list,1) 1]),reshape(spin_S, [9*size(pair_list,1) 1]),...
         reshape(oper_L, [9*size(pair_list,1) 1]),reshape(oper_S, [9*size(pair_list,1) 1]),...
         reshape(H_coeff,[9*size(pair_list,1) 1]),reshape(T_coeff,[9*size(pair_list,1) 5]),...
         reshape(irrcomp,[9*size(pair_list,1) 5]),'VariableNames',{'L','S','opL','opS','H','T','phi'});

% Kill insignificant rows
tol1=spin_system.tols.inter_cutoff; tol2=eps;
D2((abs(D2.H)<tol1)&((sum(abs(D2.phi),2)<tol1)|(sum(abs(D2.T),2)<tol2)),:)=[];

% Merge descriptors and clean up
descr=[D1; D2]; clear('D1','D2'); nterms=size(descr,1);
clear('spin_L','spin_S','oper_L','oper_S','H_coeff','T_coeff','irrcomp');
report(spin_system,[num2str(nterms) ' unique operators in the Hamiltonian descriptor.']);

% Balance the descriptor
report(spin_system,'balancing descriptor for parallel processing...');
descr=descr(randperm(size(descr,1)),:);

% Build the Hamiltonian
report(spin_system,'building the Hamiltonian...');
spmd
    
    % Localize the problem at the nodes
    partition=codistributor1d.defaultPartition(nterms);
    codistrib=codistributor1d(1,partition,[nterms 1]);
    local_terms=getLocalPart(codistributed((1:nterms)',codistrib));
    
    % Preallocate the local Hamiltonian
    H=0*unit_oper(spin_system);
    if build_aniso
        Q=cell(5,5);
        for m=1:5
            for k=1:5
                Q{k,m}=0*unit_oper(spin_system);
            end
        end
    end
    
    % Build the local Hamiltonian
    for n=local_terms'
        
        % Compute operator from the specification
        if descr.S(n)==0
            oper=operator(spin_system,descr.opL(n),{descr.L(n)},operator_type);
        else
            oper=operator(spin_system,[descr.opL(n),descr.opS(n)],{descr.L(n),descr.S(n)},operator_type);
        end
        
        % Add to relevant local arrays
        H=H+descr.H(n)*oper;
        if build_aniso
            for m=1:5
                for k=1:5
                    if abs(descr.T(n,k)*descr.phi(n,m))>spin_system.tols.inter_cutoff
                        Q{k,m}=Q{k,m}+descr.T(n,k)*descr.phi(n,m)*oper;
                    end
                end
            end
        end
        
    end
            
    % Collect the result
    H=gplus(H,1); if build_aniso, Q=gplus(Q,1); end
    
end

% Retrieve the result
H=H{1}; if build_aniso, Q=Q{1}; end

% Clean up the result
H=clean_up(spin_system,H,spin_system.tols.liouv_zero);
if build_aniso
    for m=1:5
        for k=1:5
            Q{k,m}=clean_up(spin_system,Q{k,m},spin_system.tols.liouv_zero);
        end
    end
end

% Remind the user about the anisotropic part
if ~build_aniso
    report(spin_system,'WARNING - only the isotropic part has been returned.');
end

end

% Consistency enforcement
function grumble(spin_system)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~isfield(spin_system.inter.coupling,'strength'))||...
   (~isfield(spin_system.inter.zeeman,'strength'))     
    error('assumption information is missing, run assume() before calling this function.');
end
end

% IK's favourite composer is Jeremy Soule -- much of Spinach 
% coding was done with Skyrim, Guild Wars and Oblivion sound-
% tracks in the background. It is the sign of the times that
% a worthy successor to Ennio Morricone has emerged from the
% unlikeliest of locations: videogame music.

