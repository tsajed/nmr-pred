% Performs TT-truncation (right-to-left) for a tensor train.
% Normally you should not call this subroutine directly.
%
%  TTOUT=TTCLASS_SVD(TT)
%
% TT should contain a single tensor train, i.e. TT.ntrains=1.
% TT should be orthogonalised left-to-right, eg. by TTCLASS_ORT.
% TTOUT contains the same data approximated with lower TT-ranks.
% Approximation threshold is read from TT.tolerance (in Frobenius norm).
% TTOUT is return orthogonalised right-to-left.
%
% d.savostyanov@soton.ac.uk
%
function ttout=ttclass_svd(tt)

% Read tensor ranks and dimensions
sz=tt.sizes;
rnk=tt.ranks;
d=tt.ncores;
N=tt.ntrains;

if N>1
    error('Please sum all tensor trains before calling ttclass_svd');
end

% Preallocate the result
ttout=tt;
r=rnk(:,1);

% Define the relative approximation accuracy
eps=tt.tolerance(1,1)/tt.coeff(1,1);
eps=eps/sqrt(d);

% SVD cores right-to-left
for k=d:-1:2
    C=reshape(ttout.cores{k,1}, [r(k), sz(k,1)*sz(k,2)*r(k+1)]);
    B=reshape(ttout.cores{k-1,1}, [r(k-1)*sz(k-1,1)*sz(k-1,2), r(k)]);
    
    [U,S,V]=svd(C,'econ');
    S=diag(S);
    rnew=ttclass_chop(S,eps*norm(S,'fro'));
    U=U(:,1:rnew); S=S(1:rnew); V=V(:,1:rnew);
    
    ttout.cores{k,1}=reshape(V', [rnew, sz(k,1), sz(k,2), r(k+1)]);
    ttout.cores{k-1,1}=reshape(B*U*diag(S), [r(k-1), sz(k-1,1), sz(k-1,2), rnew]);
    r(k)=rnew;
end

nrm=norm(ttout.cores{1,1}(:));
ttout.cores{1,1}=ttout.cores{1,1}/nrm;
ttout.coeff(1,1)=ttout.coeff(1,1)*nrm;
end