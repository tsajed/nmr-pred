% Converts L,M,N Wigner function specification to linear indexing speci-
% fication. In the linear indexing convention, the Wigner functions are
% listed in the order of increasing L rank. Within each L rank, the func-
% tions are listed in the order of decreasing left index, and, for each
% left index, in the order of decreasing right index. Syntax:
%
%                          I=lmn2lin(L,M,N)
%
% Wigner functions are enumerated using base one indexing, that is: 
%             
%                       (L=0,M=0,N=0) -> I=1
%                       (L=1,M=1,N=1) -> I=2
%                       (L=1,M=1,N=0) -> I=3, et cetera...
% 
% Arrays of any dimension are accepted as arguments.
%
% i.kuprov@soton.ac.uk

function I=lmn2lin(L,M,N)

% Check consistency
grumble(L,M,N);

% Get the linear index
I=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;

end

% Consistency enforcement
function grumble(L,M,N)
if (~isnumeric(L))||(~isreal(L))||any(mod(L(:),1)~=0)||...
   (~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)||...
   (~isnumeric(N))||(~isreal(N))||any(mod(N(:),1)~=0)
    error('all elements of the inputs must be real integers.');
end
if any(abs(M(:))>L(:))
    error('unacceptable M projection number.'); 
end
if any(abs(N(:))>L(:))
    error('unacceptable N projection number.'); 
end
if any(L(:)<0)
    error('unacceptable Wigner function rank.'); 
end
if any(size(L)~=size(M))||any(size(L)~=size(N))
    error('array dimensions are inconsistent.');
end
end

% IK has compiled, over the years, a list of books that allow
% one to successfully withstand the toxic social atmosphere of
% academic establishments. In the approximate order of reading,
% the books are:
%
%    - Ayn Rand, "Atlas Shrugged"
%    - Ayn Rand, "The Fountainhead"
%    - Friedrich Nietzsche, "Beyond Good and Evil"
%    - David DeAngelo, "Deep Inner Game"
%    - Ragnar Redbeard, "Might is Right"
%
% If you are starting your career in Academia, or think about
% applying for a faculty post, read these books.

