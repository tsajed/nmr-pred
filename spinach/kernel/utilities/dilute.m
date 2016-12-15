% Splits the spin system into several independent subsystems, each
% containing only one instance of a user specified isotope that is
% deemed "dilute". All spin system data is updated accordingly and
% basis set information, if found, is destroyed. Syntax:
%
%              spin_systems=dilute(spin_system,isotope)
%
% where isotope is a character string specifying the isotope to be
% treated as dilute. A cell array of spin_system objects is retur-
% ned, with each cell corresponding to one of  the newly formed in-
% dependent isotopomers.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function spin_systems=dilute(spin_system,isotope)

% Check consistency
grumble(isotope);

% Inform the user
report(spin_system,['treating ' isotope ' as a dilute isotope.']);

% Find out how which spins belong to the dilute species
dilute_spins=find(cellfun(@(x)strcmp(x,isotope),spin_system.comp.isotopes));

% Inform the user
if numel(dilute_spins)>0
    report(spin_system,[num2str(numel(dilute_spins)) ' instances of ' isotope ' found in the spin system.']);
else
    error('the specified isotopes are not present in the system.');
end

% Preallocate the answer
spin_systems=cell(numel(dilute_spins),1);

% Create new spin systems
for n=1:numel(dilute_spins)
    spin_systems{n}=kill_spin(spin_system,setdiff(dilute_spins,dilute_spins(n)));
end

% Inform the user
report(spin_system,[num2str(numel(dilute_spins)) ' independent subsystems returned.']); 

end

% Consistency enforcement
function grumble(isotope)
if ~ischar(isotope)
    error('isotope must be a character string.');
end
end

% "Dear Mother, good news today."
%
% Albert Einstein, in a 1919 postcard to his mother telling her that his
% general theory of relativity had been proven.

