% Returns the "carrier" Hamiltonian - the part of the Zeeman 
% interaction Hamiltonian that corresponds to all particles
% having the Zeeman frequency prescribed by their magnetogy-
% ric ratio and the magnet field specified by the user. This
% Hamiltonian is frequently used in rotating frame transfor-
% mations and average Hamiltonian theories. Syntax:
%
%               H=carrier(spin_system,spins)
%
% i.kuprov@soton.ac.uk

function H=carrier(spin_system,spins)

% Preallocate the answer
H=mprealloc(spin_system,1);

% Find the spins
if strcmp(spins,'all')
    spin_index=1:spin_system.comp.nspins;
else
    spin_index=find(strcmp(spins,spin_system.comp.isotopes));
end

% Compute the answer
for n=1:numel(spin_index) %#ok<*PFBNS>
    if significant(spin_system.inter.basefrqs(spin_index(n)),spin_system.tols.liouv_zero)
        H=H+spin_system.inter.basefrqs(spin_index(n))*operator(spin_system,{'Lz'},{spin_index(n)});
    end
end

% Clean up the answer
H=clean_up(spin_system,(H+H')/2,spin_system.tols.liouv_zero);

end

% Atelophobia - (n.) fear of imperfection

