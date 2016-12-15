% Converts linear indexing state specification to L,M indexing. In
% the linear indexing convention, the states are listed in the order
% of increasing L rank, and, within ranks, in the order of decreas-
% ing M projection. Syntax: [L,M]=lin2lm(I)
%
% WARNING - zero base indexing, that is: 
%             
%                           I=0 -> (L=0,M=0)
%                           I=1 -> (L=1,M=1)
%                           I=2 -> (L=1,M=0), et cetera...
% 
% Arrays of any dimension are accepted as arguments.
%
% i.kuprov@soton.ac.uk

function [L,M]=lin2lm(I)

% Check consistency
grumble(I);

% Get the ranks and projections
L=fix(sqrt(I)); M=L.^2+L-I;

% Make sure the conversion is correct
if nnz(lm2lin(L,M)~=I)>0
    error('IEEE arithmetic breakdown, please contact the developer.');
end

end

% Consistency enforcement
function grumble(I)
if (~isnumeric(I))||(~isreal(I))||any(mod(I(:),1)~=0)||any(I(:)<0)
    error('all elements of the input array must be non-negative integers.');
end
end

% Arrogance on the part of the meritorious is even more offensive
% to us than the arrogance of those without merit: for merit itself
% is offensive.
%
% Friedrich Nietzsche

