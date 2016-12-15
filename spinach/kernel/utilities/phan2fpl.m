% Projects a spatial intensity distribution into the Fokker-Planck
% space, using the spin state supplied. Syntax:
%
%                      rho=phan2fpl(phan,rho)
%
% Parameters:
%
%        phan  - phantom (the spatial distribution of the
%                amplitude of the specified spin state)
%
%        rho   - the spin state in question
%
% i.kuprov@soton.ac.uk

function rho=phan2fpl(phan,rho)

% Stretch the phantom and kron it into the spin state
rho=kron(phan(:),rho);

end

% Q: "How many members of a certain demographic group does
%     it take to perform a specified task?"
%
% A: "A finite number: one to perform the task and the rem-
%     ainder to act in a manner stereotypical of the group
%     in question."

