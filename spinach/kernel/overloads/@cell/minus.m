% Subtracts cell arrays element-by-element.
%
% i.kuprov@soton.ac.uk

function A=minus(A,B)

if iscell(A)&&iscell(B)&&all(size(A)==size(B))
    for n=1:numel(A), A{n}=A{n}-B{n}; end
else
    error('cell array sizes must match.');
end

end

% I can, therefore I am.
%
% Simone Weil

