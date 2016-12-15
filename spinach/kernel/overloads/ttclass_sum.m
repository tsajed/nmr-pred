% This subroutine packs all trains from buffer to a single train.
% Normally you should not call it directly, use SHRINK instead.
%
%  TTOUT=TTCLASS_SUM(TT)
%
% TT contains several buffered trains,
% TTOUT has one train with sum-of-ranks.
%
% d.savostyanov@soton.ac.uk
%
function ttout=ttclass_sum(tt)

% Read tensor ranks and dimensions
sz=tt.sizes;
rnk=tt.ranks;
[d,N]=size(tt.cores);

% Fast return if possible
if N==1
    ttout=tt;
    return
end

% Total rank of all summands
R=sum(rnk,2);
R(1)=1; R(d+1)=1;

% Preallocate a single tensor train for all summands
core=cell(d,1);
for k=1:d
    core{k}=zeros(R(k),sz(k,1),sz(k,2),R(k+1));
end

% Collect buffered trains in a single tensor train
rr=[zeros(d+1,1), cumsum(rnk,2)];
for n=1:N
    if d>1
        core{1,1}(1,:,:,rr(2,n)+1:rr(2,n+1))=tt.coeff(1,n)*tt.cores{1,n};
        for k=2:d-1
            core{k,1}(rr(k,n)+1:rr(k,n+1), :,:, rr(k+1,n)+1:rr(k+1,n+1))=tt.cores{k,n};
        end     
        core{d,1}(rr(d,n)+1:rr(d,n+1),:,:,1)=tt.cores{d,n};
    else
        % Special case of one-site tensor trains
        core{1,1}(1,:,:,1)=core{1,1}(1,:,:,1)+tt.coeff(1,n)*tt.cores{1,n};
    end
end

% Construct empty ttclass array and fill the fields
ttout=ttclass;
ttout.coeff=1;
ttout.cores=core;
ttout.tolerance=abs(sum(tt.tolerance));
end