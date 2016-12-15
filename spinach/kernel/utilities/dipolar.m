% Computes dipolar couplings in the presence or absence of periodic
% boundary conditions. This is an auxiliary function of Spinach ker-
% nel, direct calls are discouraged. Use xyz2dd to convert Cartesian
% coordinates to dipolar couplings.
%
% i.kuprov@soton.ac.uk

function spin_system=dipolar(spin_system)

% Report to the user
report(spin_system,['dipolar interaction distance threshold: ' num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);
report(spin_system,'running dipolar interaction network analysis...');

% Preallocate distance vector array
distvects=cell(spin_system.comp.nspins,spin_system.comp.nspins);

% Loop over chemical species
for m=1:numel(spin_system.chem.parts)
    
    % Extract the spin list
    spin_list=spin_system.chem.parts{m};
    
    % Loop over atom pairs
    for n=spin_list
        for k=spin_list
            
            % Only proceed if coordinates are specified
            if (~isempty(spin_system.inter.coordinates{n}))&&(~isempty(spin_system.inter.coordinates{k}))&&(n~=k)
                
                % Determine possible distance vectors from n to k
                if isempty(spin_system.inter.pbc)
                    
                    % Compute the distance vector for standalone system
                    dv=spin_system.inter.coordinates{k}-spin_system.inter.coordinates{n};
                    
                elseif numel(spin_system.inter.pbc)==1
                    
                    % Preallocate distance vector array
                    nvecs=2*spin_system.tols.dd_ncells+1; dv=zeros(nvecs,3);
                    
                    % Compute distance vectors with 1D periodic boundary
                    linear_index=1;
                    for p=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells
                        dv(linear_index,:)=spin_system.inter.coordinates{k}+p*spin_system.inter.pbc{1}-...
                                           spin_system.inter.coordinates{n};
                        linear_index=linear_index+1;
                    end
                    
                elseif numel(spin_system.inter.pbc)==2
                    
                    % Preallocate distance vector array
                    nvecs=(2*spin_system.tols.dd_ncells+1)^2; dv=zeros(nvecs,3);
                    
                    % Compute distance vectors with 2D periodic boundary
                    linear_index=1;
                    for p=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells
                        for q=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells
                            dv(linear_index,:)=spin_system.inter.coordinates{k}+p*spin_system.inter.pbc{1}+...
                                               q*spin_system.inter.pbc{2}-spin_system.inter.coordinates{n};
                            linear_index=linear_index+1;                
                        end
                    end
                    
                elseif numel(spin_system.inter.pbc)==3
                    
                    % Preallocate distance vector array
                    nvecs=(2*spin_system.tols.dd_ncells+1)^3; dv=zeros(nvecs,3);
                    
                    % Compute distance vectors with 3D periodic boundary
                    linear_index=1;
                    for p=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells
                        for q=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells
                            for r=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells
                                dv(linear_index,:)=spin_system.inter.coordinates{k}+p*spin_system.inter.pbc{1}+...
                                                   q*spin_system.inter.pbc{2}+r*spin_system.inter.pbc{3}-...
                                                   spin_system.inter.coordinates{n};
                                linear_index=linear_index+1;
                            end
                        end
                    end
                    
                else
                    
                    % Complain and bomb out
                    error('PBC translation vector array has invalid dimensions.');
                    
                end
                
                % Ignore distance vectors that are longer than the threshold
                dv=dv(sqrt(sum(dv.^2,2))<spin_system.tols.prox_cutoff,:);
                
                % Detect atomic collisions
                if any(sqrt(sum(dv.^2,2))<0.5)
                    
                    % Bomb out for distances closer than 0.5 Angstrom
                    error(['collision detected between spin ' num2str(n) ' and spin ' num2str(k) ' or their PBC images.']);
                    
                end
                
                % Assign the cell array
                if ~isempty(dv), distvects{n,k}=dv; end
                
            end
        end
    end
    
end

% Get proximity matrix
spin_system.inter.proxmatrix=sparse(~cellfun(@isempty,distvects));

% Find interacting atom pairs
[rows,cols,~]=find(spin_system.inter.proxmatrix);

% Report to the user
report(spin_system,['found ' num2str(numel(rows)/2) ' spin pair(s) with distance under the '...
                    'threshold of ' num2str(spin_system.tols.prox_cutoff) ' Angstrom.']);

% Loop over interacting atom pairs
for n=1:numel(rows)
    
    % Loop over PBC directions
    for k=1:size(distvects{rows(n),cols(n)},1)
        
        % Get the distance
        distance=norm(distvects{rows(n),cols(n)}(k,:),2);
        
        % Get the ort
        ort=distvects{rows(n),cols(n)}(k,:)/distance;
        
        % Compute the dipolar interaction prefactor (0.5 due to double counting)
        A=0.5*spin_system.inter.gammas(rows(n))*spin_system.inter.gammas(cols(n))*...
              spin_system.tols.hbar*spin_system.tols.mu0/(4*pi*(distance*1e-10)^3);
        
        % Compute the dipolar coupling matrix
        D=A*[1-3*ort(1)*ort(1)   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
              -3*ort(2)*ort(1)  1-3*ort(2)*ort(2)   -3*ort(2)*ort(3);
              -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)*ort(3)];
        
        % Add to the total
        spin_system.inter.coupling.matrix{rows(n),cols(n)}=spin_system.inter.coupling.matrix{rows(n),cols(n)}+D;
        
    end
    
end

end

% What worries me about religion is that it teaches people to
% be satisfied with not understanding. 
%
% Richard Dawkins

