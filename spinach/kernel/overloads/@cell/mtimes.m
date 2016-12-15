% Element-by-element multiplication of cell arrays of matrices.
%
% i.kuprov@soton.ac.uk

function C=mtimes(A,B)

if iscell(A)&&isnumeric(B)
    C=cell(size(A)); for n=1:numel(A), C{n}=A{n}*B; end
elseif isnumeric(A)&&iscell(B)
    C=cell(size(B)); for n=1:numel(B), C{n}=A*B{n}; end
else
    error('at least one argument must be numeric.');
end

end

% We keep blaming Comrade Stalin and I agree that we have good reasons. But
% I would like to ask: who wrote those four million denunciation letters?
% 
% S.D. Dovlatov

