% Directional derivatives of the matrix exponential. Implements Equation 11
% of Najfeld and Havel and our Equations 16 and 17. Syntax:
%
%                        D=dirdiff(spin_system,A,B,T,N)
%
% if A an B are matrices, computes the N-th derivative of the matrix expo-
% nential exp(A*T) in the direction B. If B is a cell array, computes the
% mixed derivative of exp(A*T) in directions B{1}, B{2}, etc.
%
% i.kuprov@soton.ac.uk
% d.goodwin@soton.ac.uk

function D=dirdiff(spin_system,A,B,T,N)

% Check consistency
grumble(A,B,T,N);

% Preallocate arrays
auxmat=cell(N,N); D=cell(1,N);

% Build auxiliary matrix
for n=1:N
    for k=1:N
        auxmat{n,k}=sparse(size(A,1),size(A,2));
    end
end
if iscell(B)
    for n=1:(N-1), auxmat{n,n+1}=B{n}; end
else
    for n=1:(N-1), auxmat{n,n+1}=B; end
end
for n=1:N, auxmat{n,n}=A; end

% Tighten up exponentiation tolerance
spin_system.tols.prop_chop=1e-14;

% Exponentiate auxiliary matrix
auxmat=propagator(spin_system,cell2mat(auxmat),T);

% Extract directional derivatives
for n=1:N
    D{n}=factorial(n-1)*auxmat(1:size(A,1),(1:size(A,2))+size(A,2)*(n-1));
end

end

% Consistency enforcement
function grumble(A,B,T,N)
if (~isnumeric(N))||(~isreal(N))||(~isscalar(N))||(N<2)||(mod(N,1)~=0)
    error('N must be a real integer greater than 1.');
end
if iscell(B)
    if numel(B)~=N-1
        error('number of B matrices must equal N-1.');
    else
        for n=1:numel(B)
            if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
                (~isnumeric(B{n}))||(size(B{n},1)~=size(B{n},2))
                error('A and B must be square matrices.');
            end
        end
    end
else
    if (~isnumeric(A))||(size(A,1)~=size(A,2))||...
        (~isnumeric(B))||(size(B,1)~=size(B,2))
        error('A and B must be square matrices.');
    end
end
if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))
    error('T must be a real scalar.');
end
end

% Disciplining yourself to do what you know is right 
% and important, although difficult, is the high road
% to pride, self-esteem, and personal satisfaction.
%
% Margaret Thatcher

