% Wipes the specified rows and columns.
%
% i.kuprov@soton.ac.uk

function spec=killcross(spec,f1idx,f2idx)

% Wipe the indices
spec(f2idx,:)=0;
spec(:,f1idx)=0;

end

% A narcissist is someone better-looking than you are.
%
% Gore Vidal

