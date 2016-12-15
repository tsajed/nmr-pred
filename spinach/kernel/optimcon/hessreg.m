% Hessian Regularisation for quasi-Newton and Newton-Raphson optimisation
% regularisations. To be used when ill-conditioned or near-singular Hessian
% matrices are expected. Can also be used to condition any square,
% symmetric matrix. 
%
% <http://spindynamics.org/wiki/index.php?title=Hessreg.m>

function [H,g,data]=hessreg(H,g,reg_param,data)

% Bootstrap if necessary
if ~exist('reg_param','var'), reg_param=bootstrap(); end

% Check consistency
grumble(H);

%initialise regularisation counter
ind=0;

switch reg_param.method
    
    case {''}
        % nothing to do
    case {'reset'}
        % simple reset hessian to the identity i.e. gradient descent iterate
        
        ind=ind+1; % increase counter
        
        %reset the hessian to identity
        H=eye(size(H));
        
    case {'CHOL'}
        % If specified, perform Cholesky regularisation
        
        if min(diag(H))>=0; sigma=norm(H,'fro');
        else sigma=norm(H,'fro')-min(diag(H)); end
        
        ind=ind+1; % increase counter
        [~,pos_def_test]=chol(H);
        while ind<reg_param.n_reg && logical(pos_def_test)
            
            [L,pos_def_test]=chol(H+sigma*eye(size(H))); iL=inv(L.');
            %output the inverse hessian
            H=iL.'*iL;
            % increase sigma and counter
            sigma=max(2*sigma,norm(H,'fro')); ind = ind+1;
            
        end
        
    case {'TRM'}
        % Trust Region Method of regularisation
        
        switch reg_param.cond_method
            
            case {''}
                
                ind=ind+1; % increase counter
                
                % eigendecomposition
                [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,H);
                
                % use the eigendecomposition to reform hessian
                H=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma);
                
            case {'iterative','scaled'}
                ind=ind+1; % increase counter
                
                reg_param.delta=1; % initialise alpha
                
                % eigendecomposition
                [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,H);
                
                if ismember(reg_param.cond_method,{'scaled'})
                    
                    % scale alpha to smallest eigenvalue
                    reg_param.delta=sqrt(abs(min(diag(hess_eigs))));
                    
                    % initial scaling of the Hessian
                    [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,H);
                end
                
                % use the eigendecomposition to reform hessian
                H=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma);
                hess_comp=H;
                
                % iteratively condition the hessian until < limit
                while cond(H,2)>reg_param.max_cond...
                        && ind<reg_param.n_reg
                    
                    % define new alpha
                    reg_param.delta=reg_param.delta*reg_param.phi;
                    
                    % eigendecomposition
                    [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,hess_comp);
                    
                    % use the eigendecomposition to reform hessian
                    H=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma);
                    
                    ind=ind+1; % increase counter
                end
        end
        
    case {'RFO'}
        % Rational Function Optimisation
        
        switch reg_param.cond_method
            
            case {''}
                
                ind=ind+1; % increase counter
                
                % eigendecomposition
                [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,H,g);
                
                % use the eigendecomposition to reform hessian
                [H,g]=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma,g);
                
            case {'iterative','scaled'}
                
                ind=ind+1; % increase counter
                
                reg_param.alpha=1; % initialise alpha
                
                % eigendecomposition
                [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,H,g);
                
                if ismember(reg_param.cond_method,{'scaled'})
                    
                    % scale alpha to smallest eigenvalue
                    reg_param.alpha=sqrt(1/abs(min(diag(hess_eigs))));
                    
                    % initial scaling of the Hessian
                    [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,H,g);
                end
                
                % use the eigendecomposition to reform hessian
                [H,g]=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma,g);
                hess_comp=H; grad_comp=g;
                
                % iteratively condition the hessian until < limit
                while cond(H,2)>reg_param.max_cond...
                        && ind<reg_param.n_reg
                    
                    % define new alpha
                    reg_param.alpha=reg_param.alpha*reg_param.phi;
                    
                    % eigendecomposition
                    [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,hess_comp,grad_comp);
                    
                    % use the eigendecomposition to reform hessian
                    [H,g]=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma,grad_comp);
                    
                    ind=ind+1; % increase counter
                end
        end
        
    otherwise
        error('unknown regularisation method')
end

% compile diagnostics data
if nargout>2 && exist('data','var')
    if ismember(reg_param.method,{'TRM','RFO'})
        [hess_eigvecs,hess_eigs]=eig(H);
        data.hess_eigs=diag(hess_eigs);
        data.hess_eigvecs=hess_eigvecs;
        data.hess_mineig=min(abs(data.hess_eigs));
        data.hess_cond=max(abs(data.hess_eigs))/data.hess_mineig;
    end
    data.n_cond=ind;
else data=[];
end

end

function [hess_eigs,hess_eigvecs,sigma]=find_eigs(reg_param,hess,grad)

switch reg_param.method
    case {'RFO'}
        
        % define the auxiliary Hessian
        hess_aug=[(reg_param.alpha^2).*hess   reg_param.alpha.*grad;
            reg_param.alpha.*grad'        0];
        
        % eigendecomposition of auxiliary Hessian
        [hess_eigvecs,hess_eigs]=eig(hess_aug);
        
        % find lambda smaller than lowest Hessian negative eigenvalue
        sigma=reg_param.def*min([0 min(reg_param.def*diag(hess_eigs))]);
        
    case {'TRM'}
        
        % eigendecompositon of Hessian
        [hess_eigvecs,hess_eigs]=eig(hess);
        
        % find lambda smaller than lowest Hessian negative eigenvalue
        sigma=reg_param.def*max([0 reg_param.delta-min(reg_param.def*diag(hess_eigs))]);
end

end

function [hess,grad]=reconstruct_hess(reg_param,hess_eigs,hess_eigvecs,sigma,grad)

% set the identity
I=eye(size(hess_eigs));

% switch between regularisation optim_param.regularisation
switch reg_param.method
    case {'RFO'}
        
        % scale Hessian to have all positive eigenvalues
        hess=hess_eigvecs*((hess_eigs-sigma.*I)/hess_eigvecs);
        
        % rescale the Hessian and gradient
        grad=hess(1:end-1,end)./(reg_param.alpha);
        
        % set Hessian to correct size
        hess=hess(1:end-1,1:end-1)./(reg_param.alpha^2);
        
    case {'TRM'}
        
        % scale Hessian to have all positive eigenvalues
        hess=hess_eigvecs*((hess_eigs+sigma.*I)/hess_eigvecs);
end

end

% Set default options if not already set
function reg_param=bootstrap()

reg_param.method='RFO';
reg_param.cond_method='iterative';
reg_param.alpha=1;
reg_param.delta=1;
reg_param.max_cond=1e4;
reg_param.phi=0.9;
reg_param.n_reg=2500;

end

% Consistency enforcement
function grumble(H)
if (~isnumeric(H))||(size(H,1)~=size(H,2))||...
        (~issymmetric(H))||(~isreal(H))
    error('The Hessian must be a symmmetric, real, square, matrix.');
end
end

% A youth who had begun to read geometry with Euclid, when he had learnt
% the first proposition, inquired, "What do I get by learning these
% things?". So Euclid called a slave and said "Give him three pence, since
% he must make a gain out of what he learns."
%
% Joannes Stobaeus
