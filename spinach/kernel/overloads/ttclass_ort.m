% Performs TT-orthogonalisation for all trains from buffer.
% DIR=+1 {default} gives you left-to-right orthogonality,
% DIR=-1 right-to-left
% Normally you should not call this subroutine directly.
%
%  TT=TTCLASS_ORT(TT,DIR)
%  [TT,LOGNRM]=TTCLASS_ORT(TT,DIR)
%
% On input, TT contains several buffered trains,
% On output, TT has all of them orthogonalised left-to-right.
% If LOGNRM is present, all buffered trains are also normalized,  
% and natural logarithms of their norms returned in the vector LOGNRM.
% Use this option if the tensor norm is likely to exceed realmax()=1.7977e+308. 
%
% d.savostyanov@soton.ac.uk
%
function [tt,lognrm]=ttclass_ort(tt,dir)

% Read tensor ranks and dimensions
sz=tt.sizes;
rnk=tt.ranks;
d=tt.ncores;
N=tt.ntrains;

% Set the default dir if it is not present
if ~exist('dir','var')
    dir=+1;
end
if (nargout>1)
    lognrm=log(tt.coeff);
    tt.coeff=ones(1,N);
end

switch dir
    case +1
        for n=1:N
            % TT-othogonalise each tensor train left-to-right
            r=rnk(:,n);
            for k=1:d-1
                B=reshape(tt.cores{k,n}, [r(k)*sz(k,1)*sz(k,2), r(k+1)]);
                C=reshape(tt.cores{k+1,n}, [r(k+1), sz(k+1,1)*sz(k+1,2)*r(k+2)]);
                
                [Q,R]=qr(B,0);
                nrm=norm(R);
                
                % Take care of the norm
                if (nrm>0)
                    R=R/nrm;
                    if (nargout>1)
                        lognrm(1,n)=lognrm(1,n)+log(nrm);
                    else
                        tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                    end
                end
                C=R*C;
                rnew=size(Q,2);
                
                tt.cores{k,n}=reshape(Q, [r(k), sz(k,1), sz(k,2), rnew]);
                tt.cores{k+1,n}=reshape(C, [rnew, sz(k+1,1), sz(k+1,2), r(k+2)]);
                r(k+1)=rnew;
            end
            
            % Take care of the last norm
            nrm=norm(tt.cores{d,n}(:));
            if (nrm>0)
                tt.cores{d,n}=tt.cores{d,n}/nrm;
                if (nargout>1)
                    lognrm(1,n)=lognrm(1,n)+log(nrm);
                else
                    tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                end
            end
        end
    case -1
        for n=1:tt.ntrains
            % TT-othogonalise each tensor train lright-to-left
            r=rnk(:,n);
            for k=d:-1:2
                B=reshape(tt.cores{k,n}, [r(k), sz(k,1)*sz(k,2)*r(k+1)]);
                C=reshape(tt.cores{k-1,n}, [r(k-1)*sz(k-1,1)*sz(k-1,2), r(k)]);
                
                [Q,R]=qr(B.',0);
                nrm=norm(R);
                
                % Take care of the norm
                if (nrm>0)
                    R=R/nrm;
                    if (nargout>1)
                        lognrm(1,n)=lognrm(1,n)+log(nrm);
                    else
                        tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                    end
                end
                C=C*R.';
                rnew=size(Q,2);
                
                tt.cores{k,n}=reshape(Q.', [rnew, sz(k,1), sz(k,2), r(k+1)]);
                tt.cores{k-1,n}=reshape(C, [r(k-1), sz(k-1,1), sz(k-1,2), rnew]);
                r(k)=rnew;
            end
            
            % Take care of the last norm
            nrm=norm(tt.cores{1,n}(:));
            if (nrm>0)
                tt.cores{1,n}=tt.cores{1,n}/nrm;
                if (nargout>1)
                    lognrm(1,n)=lognrm(1,n)+log(nrm);
                else
                    tt.coeff(1,n)=tt.coeff(1,n)*nrm;
                end
            end
        end
    otherwise
        error('unrecognized parameter dir.');
end
end