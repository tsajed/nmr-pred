% Cluster Correlation Expansion function - partitions the spin_system
% object into individual subsystems prescribed by the CCE clusterisa-
% tion model [DOI reference here]. Syntax:
%
%       [subsystems,indices]=cce(spin_system,system,bath,order)
%
% Parameters:
%
%     system   - a row vector giving the numbers of the
%                spins belonging to the "system"
%
%     bath     - a row vector giving the numbers of the
%                spins belonging to the "bath"
%
%     order    - CCE expansion order, the number of bath
%                spins in each subsystem
%
% Outputs:
%
%     subsystems - a cell array of spin system objects
%                  for the individual subsystems
%
%     indices    - a cell array of index vectors indi-
%                  cating which bath spins ended up in
%                  each of the generated subsystems.
%                  The spins are listed in the order
%                  of increasing index.
%
% i.kuprov@soton.ac.uk
% lipiongyang@csrc.ac.cn

function [subsystems,indices]=cce(spin_system,system,bath,order)

% Start the index array
indices={};

% Loop over bath cluster sizes
for n=1:order
    
    % Generate bath spin picks
    partial_index=combnk(bath,n);
    
    % Sort in ascending order
    partial_index=sort(partial_index,2,'ascend');
    
    % Remove repetitions
    partial_index=unique(partial_index,'rows');
    
    % Convert into a cell array
    partial_index=mat2cell(partial_index,ones(size(partial_index,1),1),n);
    
    % Add to the total
    indices=[indices; partial_index];  %#ok<AGROW>

end

% Build the corresponding subsystems
subsystems=cell(numel(indices,1));
parfor n=1:numel(indices)
    
    % Remove the spins that are not included in each cluster
    subsystems{n}=kill_spin(spin_system,setdiff(1:spin_system.comp.nspins,...
                                                [system indices{n}]));

end
    
end

