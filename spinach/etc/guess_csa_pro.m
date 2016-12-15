% Chemical shieft anisotropy estimates for peptide bond heteroatoms.
%
% <http://spindynamics.org/wiki/index.php?title=Guess_csa_pro.m>

function CSAs=guess_csa_pro(aa_nums,pdb_ids,coords)

% Check consistency
grumble(aa_nums,pdb_ids,coords);

% Preallocate CSA array
CSAs=cell(numel(pdb_ids),1);

% Number the atoms
numbers=1:numel(coords);

% Loop over amino acids
for n=2:(max(aa_nums)-1)
    
    % Assign amide nitrogen CSAs
    if ismember('C',pdb_ids(aa_nums==n))&&...
       ismember('N',pdb_ids(aa_nums==(n+1)))&&...
       ismember('H',pdb_ids(aa_nums==(n+1)))
    
        % Get C coordinates
        local_coords=coords(aa_nums==n);
        C=local_coords{strcmp('C',pdb_ids(aa_nums==n))}; C=C(:);
   
        % Get N coordinates
        local_coords=coords(aa_nums==(n+1));
        N=local_coords{strcmp('N',pdb_ids(aa_nums==(n+1)))}; N=N(:);
   
        % Get H coordinates
        local_coords=coords(aa_nums==(n+1));
        H=local_coords{strcmp('H',pdb_ids(aa_nums==(n+1)))}; H=H(:);
        
        % Get the primary directions
        N_CO_vec=C-N; N_H_vec=H-N;
        
        % Double-check the distances
        if (norm(N_CO_vec)>2.0)||(norm(N_H_vec)>2.0)
            error('Amino acid numbering is not sequential.');
        end
        
        % Make ZZ eigenvector collinear with N-CO bond
        zz_eigvec=N_CO_vec;
        zz_eigvec=zz_eigvec/norm(zz_eigvec);
        
        % Make YY eigenvector perpendicular to the peptide plane
        yy_eigvec=cross(N_CO_vec,N_H_vec);
        yy_eigvec=yy_eigvec/norm(yy_eigvec);
        
        % Make XX eigenvector perpendicular to the other two
        xx_eigvec=cross(yy_eigvec,zz_eigvec);
        xx_eigvec=xx_eigvec/norm(xx_eigvec);
        
        % Build the eigenvalue matrix
        D=diag([-125 45 80]);
        
        % Build the eigenvector matrix
        V=[xx_eigvec yy_eigvec zz_eigvec];
        
        % Identify the nitrogen
        local_numbers=numbers(aa_nums==(n+1));
        nitrogen_number=local_numbers(strcmp('N',pdb_ids(aa_nums==(n+1))));
        
        % Compose the tensor
        CSAs{nitrogen_number}=V*D*V';
        
        % Report to the user
        disp(['Amide nitrogen CSA for residue ' num2str(n) ' guessed from local geometry.']);
    
    end
    
    % Assign C=O carbon CSAs
    if ismember('C',pdb_ids(aa_nums==n))&&...
       ismember('N',pdb_ids(aa_nums==(n+1)))&&...
       ismember('CA',pdb_ids(aa_nums==n))
    
        % Get C coordinates
        local_coords=coords(aa_nums==n);
        C=local_coords{strcmp('C',pdb_ids(aa_nums==n))}; C=C(:);
   
        % Get N coordinates
        local_coords=coords(aa_nums==(n+1));
        N=local_coords{strcmp('N',pdb_ids(aa_nums==(n+1)))}; N=N(:);
   
        % Get CA coordinates
        local_coords=coords(aa_nums==n);
        CA=local_coords{strcmp('CA',pdb_ids(aa_nums==n))}; CA=CA(:);
        
        % Get the primary directions
        N_C_vec=C-N; C_CA_vec=CA-C;
        
        % Double-check the distances
        if (norm(N_C_vec)>2.0)||(norm(C_CA_vec)>2.0)
            error('Amino acid numbering is not sequential.');
        end
        
        % Make XX eigenvector collinear with C-CA bond
        xx_eigvec=C_CA_vec;
        xx_eigvec=xx_eigvec/norm(xx_eigvec);
        
        % Make ZZ eigenvector perpendicular to the >C=O plane
        zz_eigvec=cross(C_CA_vec,N_C_vec);
        zz_eigvec=zz_eigvec/norm(zz_eigvec);
        
        % Make YY eigenvector perpendicular to the other two
        yy_eigvec=cross(xx_eigvec,zz_eigvec);
        yy_eigvec=yy_eigvec/norm(yy_eigvec);
        
        % Build the eigenvalue matrix
        D=diag([70  5  -75]);
        
        % Build the eigenvector matrix
        V=[xx_eigvec yy_eigvec zz_eigvec];
        
        % Identify the carbon
        local_numbers=numbers(aa_nums==n);
        carbon_number=local_numbers(strcmp('C',pdb_ids(aa_nums==n)));
        
        % Compose the tensor
        CSAs{carbon_number}=V*D*V';
        
        % Report to the user
        disp(['Carboxyl carbon CSA for residue ' num2str(n) ' guessed from local geometry.']);
    
    end

end

end

% Consistency enforcement
function grumble(aa_nums,pdb_ids,coords)
if ~isnumeric(aa_nums)
    error('aa_nums must be a vector of integers.'); 
end
if ~iscell(pdb_ids)
    error('pdb_ids must be a cell array of character strings.');
end
if ~iscell(coords)
    error('coords must be a cell array of 3-vectors.');
end
if (numel(aa_nums)~=numel(pdb_ids))||(numel(pdb_ids)~=numel(coords))
    error('the input parameters must have the same number of elements.');
end
end

% "Where's the PCM solvent?"
%
% A 3rd year project student,
% rummaging through IK's sol-
% vent cabinet.

