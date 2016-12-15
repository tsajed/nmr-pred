% function to prepare a Hessian matrix to be definite and well-conditioned
%
% <http://spindynamics.org/wiki/index.php?title=Hessprep.m>

function [grad,hess,data]=hessprep(grad,hess,optim,data)

% force real and symmetric hessian
if ~isreal(hess), hess=real(hess); end
if ~issymmetric(hess), hess=(hess+hess')./2; end

% test for definite Hessian
if (~isempty(optim.reg.method)) || (~isempty(optim.reg.cond_method))
    
    if ismember(optim.algo.extremum,{'minimum'}) % positive definite test
        
        % check if positive definite
        [~,def_test]=chol(hess);
        
    elseif ismember(optim.algo.extremum,{'maximum'}) % negative definite test
        
        % find largest and smallest eigenvalues
        [~,eig_bnd,conv_test]=eigs(hess,2,'BE',struct('isreal',true,'issymmetric',true));
        
        if ~(logical(conv_test)), def_test=1; % eigenvalues didn't converge
        elseif all(diag(eig_bnd))>0, def_test=0; % positive definite Hessian
        else def_test=1; % all else, regularise the Hessian
        end
    end
    
    % if definite hessian, check the condition number
    if ~(logical(def_test)), hess_cond=cond(hess); end
end

if (~isempty(optim.reg.method) && (logical(def_test)) ) ||...
        (~isempty(optim.reg.cond_method) && optim.reg.max_cond<hess_cond)
    
    % regularise or/and condition the hessian
    [hess,grad,data]=hessreg(hess,grad,optim.reg,data);
    
elseif optim.reg.track_eig 
    
    % explicit tracking of Hessian eigenvalues
    [data.hess_eigvecs,data.hess_eigs]=eig(hess);
    data.hess_mineig=min(diag(data.hess_eigs));
    data.hess_cond=max(diag(data.hess_eigs))/data.hess_mineig;
    data.n_cond=0;
end

end