% Penalty terms for the Optimal Control module. Returns the penalty
% function and its gradient for the waveform, which should be sup-
% plied as row vector or a horizontal stack thereof. 
%
% <http://spindynamics.org/wiki/index.php?title=Penalty.m>

function [pen_term,pen_grad,pen_hess]=penalty(wf,type,fb,cb)

% Check consistency
if ismember(type,{'SNS','SNC'})
    grumble(wf,type,fb,cb);
else
    grumble(wf,type);
end

% Preallocate the results
if nargout>0, pen_term=0; end
if nargout>1, pen_grad=zeros(size(wf)); end
if nargout>2, pen_hess=zeros(numel(wf)); end

% Decide the penalty type
switch type
    
    case 'NS'
        
        % Compute the penalty
        if nargout>0
            pen_term=sum(sum(wf.^2));
            pen_term=pen_term/size(wf,2);
        end
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*wf;
            pen_grad=pen_grad/size(wf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*speye(numel(wf));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'DNS'
        
        % Differentiation matrix
        D=fdmat(size(wf,2),5,1);
        
        % Compute the penalty
        if nargout>0
            dwf=wf*D';
            pen_term=sum(sum(dwf.^2));
            pen_term=pen_term/size(dwf,2);
        end
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*dwf*D;
            pen_grad=pen_grad/size(dwf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*kron(D'*D,speye(size(wf,1)));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'SNS'
        
        % Build ceiling hole inventory
        ch_map=(wf>cb); ch_actual=wf.*ch_map; ch_wanted=cb.*ch_map;
        
        % Build floor hole inventory
        fh_map=(wf<fb); fh_actual=wf.*fh_map; fh_wanted=fb.*fh_map;
        
        % Compute the penalty
        pen_term=pen_term+sum(sum((ch_actual-ch_wanted).^2))+...
            sum(sum((fh_actual-fh_wanted).^2));
        pen_term=pen_term/size(wf,2);
        
        % Compute the gradient
        if nargout>1
            pen_grad=2*(ch_actual-ch_wanted)+2*(fh_actual-fh_wanted);
            pen_grad=pen_grad/size(wf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=2*ch_map/size(wf,2)+2*fh_map/size(wf,2);
            pen_hess=diag(pen_hess(:));
            pen_hess=pen_hess/size(wf,2);
        end
        
    case 'SNC'
        
        % Build ceiling hole inventory
        ch_map=(wf>cb); ch_actual=wf.*ch_map; ch_wanted=cb.*ch_map;
        
        % Build floor hole inventory
        fh_map=(wf<fb); fh_actual=wf.*fh_map; fh_wanted=fb.*fh_map;
        
        % Compute the penalty
        pen_term=pen_term+(1/6)*sum(sum((abs(ch_actual-ch_wanted)).^3))+...
            (1/6)*sum(sum((abs(fh_actual-fh_wanted)).^3));
        pen_term=pen_term/size(wf,2);
        
        % Compute the gradient
        if nargout>1
            pen_grad=(1/2)*(ch_actual-ch_wanted).^2+(1/2)*(fh_actual-fh_wanted).^2;
            pen_grad=pen_grad/size(wf,2);
        end
        
        % Compute the Hessian
        if nargout>2
            pen_hess=(ch_actual-ch_wanted)+(fh_actual-fh_wanted);
            pen_hess=diag(pen_hess(:));
            pen_hess=pen_hess/size(wf,2);
        end
        
    otherwise
        
        error('unknown penalty function type.');
        
end

end

% Consistency enforcement
function grumble(wf,type,fb,cb)
if ~isnumeric(wf)||(~isreal(wf))
    error('waveform must be a real numeric array.');
end
if ~ischar(type)
    error('type must be a character string.');
end
if ismember(type,{'SNS','SNC'})
    if ~isnumeric(fb)||(~isreal(fb))
        error('floor_bound must be a real numeric array.');
    end
    if ~isnumeric(cb)||(~isreal(cb))
        error('ceiling_bound must be a real numeric array.');
    end
    if any(size(wf)~=size(fb))
        error('floor bound array must have the same dimension as the waveform array.');
    end
    if any(size(wf)~=size(cb))
        error('ceiling bound array must have the same dimension as the waveform array.');
    end
    if any(any(cb<=fb))
        error('the ceiling bound should be above the floor bound.');
    end
end
end

% I would never die for my beliefs because I might be wrong.
%
% Bertrand Russell

