% Designates specific states as "dark" -- any population reaching
% them would end up being summed up and stored in them forever in
% a frozen state. Syntax:
%
%                L=absorb(spin_system,L,dark_states)
%
% where L is the Liouvillian and dark_states contains the numbers
% of the states that should be switched into the absorption mode.
%
% Note: absorption functionality is only available in the sphten-
%       liouv formalism.
%
% i.kuprov@soton.ac.uk

function L=absorb(spin_system,L,dark_states)

% Check consistency
grumble(spin_system,L,dark_states);

% Zero columns corresponding to dark states
L(:,dark_states)=0;

end

% Consistency enforcement
function grumble(spin_system,L,dark_states)
if ~strcmp(spin_system.bas.formalims,'sphten-liouv')
    error('this function is only applicable to sphten-liouv formalism.');
end
if (~isnumeric(L))||(size(L,1)~=size(L,2))
    error('L must be a square matrix.');
end
if any(size(L)~=size(spin_system.bas.basis,1))
    error('Dimension of the Liouvillian must match the dimension of the basis set.');
end
if (~isnumeric(dark_states))||(~isvector(dark_states))||(~isreal(dark_states))||...
   any(mod(dark_states,1)~=0)||any(dark_states<=0)
    error('dark_states must be a vector of positive integers.');
end
if any(dark_states>size(spin_system.bas.basis,1))
    error('an element of the dark_states vector exceeds the state space dimension.');
end
end

% Acting in a bad play is like spitting into eternity.
%
% Faina Ranevskaya

