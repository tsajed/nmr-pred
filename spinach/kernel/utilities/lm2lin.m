% Converts L,M spin state specification to linear indexing specification. 
% In the linear indexing convention, the states are listed in the order of
% inc reasing L rank, and, within ranks, in the order of decreasing M pro-
% jection. Syntax: I=lm2lin(L,M)
%
% WARNING - zero base indexing, that is: 
%             
%                           (L=0,M=0) -> I=0
%                           (L=1,M=1) -> I=1
%                           (L=1,M=0) -> I=2, et cetera...
%
% Arrays of any dimension are accepted as arguments.
%
% i.kuprov@soton.ac.uk

function I=lm2lin(L,M)

% Check consistency
grumble(L,M);

% Get the linear index
I=L.^2+L-M;

end

% Consistency enforcement
function grumble(L,M)
if (~isnumeric(L))||(~isreal(L))||any(mod(L(:),1)~=0)||...
   (~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)
    error('all elements of the inputs must be real integers.');
end
if any(abs(M(:))>L(:))
    error('unacceptable projection number.');
end
if any(L(:)<0)
    error('unacceptable total angular momentum.');
end
if any(size(L)~=size(M))
    error('array dimensions are inconsistent.');
end
end

% A casual stroll through the lunatic asylum shows that faith does not
% prove anything.
%
% Friedrich Nietzsche

