% Generates a state vector from its spherical tensor expansion produced
% by zeeman2sphten() function. Syntax:
%
%              rho=sphten2state(spin_system,stexp,spin_num)
%
% where stexp is a cell array with the first element of each row giving
% the operator name and the second element being the corresponding sphe-
% rical tensor expansion coefficient. The last argument is the number of
% the spin to which the expansion refers.
%
% lipingyang87@gmail.com
% i.kuprov@soton.ac.uk

function rho=sphten2state(spin_system,stexp,spin_num)

% Validate the input
grumble(spin_system,stexp,spin_num);

% Get the state vector
rho=sparse(0);
for n=1:size(stexp,1)
    
    % Get the state    
    vn=state(spin_system,stexp(n,1),{spin_num});
    
    % Normalise the state
    vn=vn/sqrt(vn(:)'*vn(:));
    
    % Add to the total
    rho=rho+stexp{n,2}*vn;
    
end

end

% Input validation function
function grumble(spin_system,stexp,spin_num)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if ~iscell(stexp)
    error('stexp must be a cell array.');
end
if (~isnumeric(spin_num))||(~isreal(spin_num))||(mod(spin_num,1)~=0)||(spin_num<1)
    error('spin_num must be a positive integer.');
end
end

