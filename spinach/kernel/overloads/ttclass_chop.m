% This subroutine defines a rank of truncated SVD decomposition 
% Normal users whould not call it directly.
%
%  R=TTCLASS_CHOP(S,TOLERANCE)
%
% Here S is a vector od singular values for a matrix,
%      TOLERANCE is absolute truncation threshold (Frobenius norm).
% Returns the optimal rank of low-rank approximation.
%
% d.savostyanov@soton.ac.uk
%
function r=ttclass_chop(s,tolerance)

x=cumsum(s(end:-1:1).^2);
k=find(x>=tolerance^2,1);
if isempty(k)
    r=0;
else
    r=numel(s)-k+1;
end

end