% Generates a Lindblad superoperator from user-specified left-side and
% right-side product superoperators and calibrates it using the experi-
% mental relaxation rate of a user-specified state. Syntax:
%
%              R=lindbladian(A_left,A_right,rho,rlx_rate)
%
% where A_left is the left side product superoperator of the interacti-
% on that is being modulated, A_right is the right side product super-
% operator of the same interaction, rho is the state vector whose rela-
% xation rate is known from the experiment and rlx_rate is that rate.
%
% i.kuprov@soton.ac.uk

function R=lindbladian(A_left,A_right,rho,rlx_rate)

% Check consistency
grumble(A_left,A_right,rho,rlx_rate);

% Generate a Lindbladian
R=A_left*A_right'-(A_left'*A_left+A_right*A_right')/2;

% Calibrate the Lindbladian
rho=rho/norm(rho,2); R=-rlx_rate*R/(rho'*R*rho);

end

% Consistency enforcement
function grumble(A_left,A_right,rho,rlx_rate)
if (~isnumeric(A_left))||(~isnumeric(A_right))||...
   (~isnumeric(rho))||(~isnumeric(rlx_rate))
    error('all inputs must be numeric.');
end
if (~isnumeric(rlx_rate))||(~isscalar(rlx_rate))||...
   (~isreal(rlx_rate))||(nsteps<0)
    error('rlx_rate must be a non-negative real number.');
end
end

% Morality, it could be argued, represents the way that people
% would like the world to work, wheareas economics represents
% how it actually does work.
%
% Steven D. Levitt, "Freakonomics"

