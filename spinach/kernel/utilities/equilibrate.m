% Equilibrates a linear chemical kinetics system for
% a user-specified time. Syntax:
%
%               c=equilibrate(K,c0,t)
%
% i.kuprov@soton.ac.uk

function c=equilibrate(K,c0,t)

% Run the chemical kinetics for a specified time
c=expm(K*t)*c0;

end

