% Computes matrix exponential integrals of the following general type:
%
%             Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]
%
% Matrix A must be Hermitian. For further information see the paper by
% Charles van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
%
%                     R=expmint(spin_system,A,B,C,T)
% 
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk

function R=expmint(spin_system,A,B,C,T)

% Check consistency
grumble(A,B,C,T);

% Build auxiliary matrix
auxmat=[-A, 1i*B; 0*A, -C];

% Compute the integral
auxmat=propagator(spin_system,auxmat,T);
R=auxmat(1:(end/2),1:(end/2))'*...
  auxmat(1:(end/2),(end/2+1):end);

% Clean up the result
R=clean_up(spin_system,R,spin_system.tols.liouv_zero);

end

% Consistency enforcement
function grumble(A,B,C,T)
if (~isnumeric(A))||(~isnumeric(B))||(~isnumeric(C))||...
   (~ismatrix(A))||(~ismatrix(B))||(~ismatrix(C))
    error('A, B and C arguments must be matrices.');
end
if (~all(size(A)==size(B)))||(~all(size(B)==size(C)))
    error('A, B and C matrices must have the same dimension.');
end
if ~ishermitian(A)
    error('A matrix must be Hermitian.');
end

if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))
    error('T must be a real scalar.');
end
end

% LADY NANCY ASTOR: "If you were my husband, Winston, I'd put poison in your tea."
% WINSTON CHURCHILL: "If I were your husband, Nancy, I'd drink it."

