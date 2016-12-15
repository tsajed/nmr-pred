% Converts D and E zero-field splitting parameters to a diagonal matrix.
%
% i.kuprov@soton.ac.uk

function M=zfs2mat(D,E)

% Compute the matrix
M=[-D/3+E, 0, 0; 0, -D/3-E, 0; 0, 0, 2*D/3];

end

% To watch the courageous Afghan freedom fighters battle modern
% arsenals with simple hand-held weapons is an inspiration to 
% those who love freedom. Their courage teaches us a great lesson - 
% that there are things in this world worth defending. To the Afghan
% people, I say on behalf of all Americans that we admire your 
% heroism, your devotion to freedom, and your relentless struggle 
% against your oppressors. 
%
% Ronald Reagan, 21 March 1983

