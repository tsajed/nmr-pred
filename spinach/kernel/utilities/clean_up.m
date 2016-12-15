% Array clean-up utility. Drops non-zero elements with magnitude below the
% user-specified tolerance and converts between sparse and dense represen-
% tations depending on the density of non-zeros in the array. Syntax:
%
%                 A=clean_up(spin_system,A,nonzero_tol)
%
% i.kuprov@soton.ac.uk
% alex_nevzorov@ncsu.edu

function A=clean_up(spin_system,A,nonzero_tol)

% Check consistency
grumble(A,nonzero_tol);

% Check if clean-up is allowed
if ~ismember('clean-up',spin_system.sys.disable)
    
    % Clean up the matrix
    A=nonzero_tol*round((1/nonzero_tol)*A);
    
    % A small matrix should always be full
    if issparse(A)&&any(size(A)<spin_system.tols.small_matrix), A=full(A); end
    
    % A big sparse matrix with too many non-zeros should be full
    if issparse(A)&&(nnz(A)/numel(A)>spin_system.tols.dense_matrix), A=full(A); end
    
    % A big full matrix with too few non-zeros should be sparse
    if (~issparse(A))&&(nnz(A)/numel(A)<spin_system.tols.dense_matrix)&&...
                       (all(size(A)>spin_system.tols.small_matrix))
        A=sparse(A);
    end
    
end

end

% Consistency enforcement
function grumble(A,nonzero_tol) %#ok<INUSL>
if ~isnumeric(nonzero_tol)
    error('the tolerance parameter must be numeric.');
end
if (~isreal(nonzero_tol))||(numel(nonzero_tol)~=1)||(nonzero_tol<=0)
    error('nonzero_tol parameter must be a positive real number.');
end
end

% The most dangerous man to any government is the man who is able to think
% things out for himself, without regard to the prevailing superstitions
% and taboos.
%
% H.L. Mencken

