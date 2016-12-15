% Krotov type algorithms for phase-sensitive state-to-state transfer prob-
% lems.
%
% <http://spindynamics.org/wiki/index.php?title=Fminkrotov.m>

function [x,fx,grad,data]=fminkrotov(cost_function,x_0,optim,hess_init,cost_fun_vars)

% Initialize the structure of the initial guess
x = x_0(:); data.x_shape=size(x_0); data.n_vars = numel(x_0);

% Always Bootstrap and set optimisation tolerances
% - assume user isn't optimisation expert
optim.cost_function=cost_function;
optim=optim_tols(optim,data.n_vars);

% add cost function variables to the data structure (empty if not provided)
if ~exist('cost_fun_vars','var'), cost_fun_vars=[]; end
data.cost_fun_vars=cost_fun_vars;

% penalty terms to consider (default=0), affects output display only.
data.npen=optim.algo.npen;

% initialise counters and timers
data.iter=0; data.n_fx=0; data.n_grad=0; data.n_cond=0; data.firstorderopt=[];
data.timeTotal=tic; data.time_fx=0; data.time_gfx=0; data.iter_timer=tic;

% print optimisation header
iterative_reporting([],'head',[],optim,data);

% initialise exitflag, functional and its derivatives
exitflag=[];

% set initial Hessian approximation
if ~exist('hess_init','var'), optim.hess=eye(data.n_vars);
else optim.hess=hess_init; end

if ismember(optim.algo.method,{'krotov-bfgs','krotov-sr1'})
    [data,x,vec_fwd,vec_bwd,fx,grad]=gradient_function(cost_function,x,[],[],data,optim.hess);
    optim.hess=real(optim.hess+optim.hess')./2;
    [grad,hess,exitflag,data]=hessian_prep(grad,optim.hess,optim.reg,data);
    optim.hess=real(hess+hess')./2;
else
    [data,x,vec_fwd,vec_bwd,fx]=gradient_function(cost_function,x,[],[],data,optim.hess);
end

% Show the current iteration information
data.itertime=toc(data.iter_timer); data.iter_timer=tic;
% report diagnostics to user
iterative_reporting(fx,'iter',exitflag,optim,data);

% Enter the main interation loop
while (true)
    
    % Keep the variables for next iteration
    store.x=x(:); store.fx=fx;
    
    % Update number of iterations
    if (isempty(exitflag)), data.iter=data.iter+1; end
    
    if ismember(optim.algo.method,{'krotov-bfgs','krotov-sr1'})
        store.grad=grad(:);
        [data,x,vec_fwd,vec_bwd,fx,grad]=gradient_function(cost_function,x,vec_fwd,vec_bwd,data,optim.hess);
    else
        [data,x,vec_fwd,vec_bwd,fx]=gradient_function(cost_function,x,vec_fwd,vec_bwd,data,optim.hess);
    end
    
    % check optimality conditions
    if (data.iter               >   optim.tols.max_iter), exitflag=0; end
    if (optim.algo.max_min*fx   <   optim.algo.max_min*store.fx), exitflag=1; end
    if (optim.algo.max_min*optim.tols.fx  < optim.algo.max_min*fx), exitflag=4; end
    if (optim.tols.x            >   max(abs(store.x(:)-x(:)))), exitflag=2; end
    
    % Check if exitflag is set from 3 previous checks
    if(~isempty(exitflag)), x=store.x; fx=store.fx; data.iter=data.iter-1; break, end;
    
    % Show the current iteration information
    data.itertime=toc(data.iter_timer); data.iter_timer=tic;
    
    % report diagnostics to user
    iterative_reporting(fx,'iter',exitflag,optim,data);
    
    % Call output function if specified
    if(output_function(optim.sys.outputfun,'iter')), exitflag=-1; end
    
    % Compute the gradient for bfgs or sr1
    if ismember(optim.algo.method,{'krotov-bfgs','ktotov-sr1'})
        
        [hess,store]=hessian_update(x(:),grad,optim.hess,store,optim.algo);
        optim.hess=real(hess+hess')./2;
        [grad,hess,data]=hessprep(grad,optim.hess,optim,data);
        optim.hess=real(hess+hess')./2;
    end
end

% Call output function if specified
if(output_function(optim.sys.outputfun,'done')), exitflag=-1; end

% Set outputs
x=reshape(x,data.x_shape);
if exist('grad','var')
    grad=reshape(grad,data.x_shape);
else
    grad=[];
end
data=iterative_reporting(fx,'term',exitflag,optim,data);

    function exit=output_function(outputfun,where)
        exit=false; % initialise
        if(~isempty(outputfun)) % data structure to output function
            % <add any variable to export here>
            % if exit is returned as true, optimisation will finish
            exit=feval(outputfun,reshape(x,data.x_shape),data,where);
        end
    end

end

function [data,x,traj_fwd,traj_bwd,fx,grad]=gradient_function(obj_fun_handle,x,vec_fwd,vec_bwd,data,hess)

% Call the cost function with gradient and Hessian if nessesary
timem=tic;
if ( nargout ==5 )
    if ~isempty(data.cost_fun_vars)
        [x,traj_fwd,traj_bwd,fx,data.cost_fun_vars]=...
            feval(obj_fun_handle,reshape(x,data.x_shape),vec_fwd,vec_bwd,data.cost_fun_vars);
    else
        [x,traj_fwd,traj_bwd,fx]=feval(obj_fun_handle,reshape(x,data.x_shape),vec_fwd,vec_bwd);
    end
    data.time_fx=data.time_fx+toc(timem);
    data.n_fx=data.n_fx+1;
elseif ( nargout ==6 )
    if ~isempty(data.cost_fun_vars)
        [x,traj_fwd,traj_bwd,fx,grad,data.cost_fun_vars]=...
            feval(obj_fun_handle,reshape(x,data.x_shape),vec_fwd,vec_bwd,hess,data.cost_fun_vars);
    else
        [x,traj_fwd,traj_bwd,fx,grad]=feval(obj_fun_handle,reshape(x,data.x_shape),vec_fwd,vec_bwd,hess);
    end
    data.time_gfx=data.time_gfx+toc(timem);
    data.n_fx=data.n_fx+1; data.n_grad=data.n_grad+1;
    grad=grad(:);
end
end

function [grad,hess,exitflag,data]=hessian_prep(grad,hess,reg,data)

% tests for positive definite and initialise exitflag
[~,pos_def_test]=chol(hess); exitflag=[];
if ~isreal(hess), exitflag =-5; % test real hessian
elseif ~issymmetric(hess), exitflag=-4; % test for symmetric hessian
elseif (~isempty(reg.method) && (logical(pos_def_test)) ) ||...
        (~isempty(reg.cond_method) && reg.max_cond<cond(hess))
    
    % regularise or/and condition the hessian
    [hess,grad,data]=hessreg(hess,grad,reg,data);
elseif reg.track_eig
    % explicit tracking of Hessian eigenvalues (for debugging)
    [data.hess_eigvecs,data.hess_eigs]=eig(hess);
    data.hess_mineig=min(diag(data.hess_eigs));
    data.hess_cond=max(diag(data.hess_eigs))/data.hess_mineig;
    data.n_cond=0;
end

end

function [hess,store] = hessian_update(x,grad,hess,store,algo)

% hessian update for quasi-Newton methods and hybrid-quasi-Newton methods
hess=quasinewton(hess, x, store.x, -algo.max_min*grad, -algo.max_min*store.grad, false, algo.method(8:end));

end

function data=iterative_reporting(fx,report_case,exitflag,optim,data)

switch report_case
    case {'head'}
        % Report start of optimisation to user
        switch(optim.algo.method)
            case 'krotov', optim_report(optim.sys.output,'start krotov optimisation',1);
            case 'krotov-sr1', optim_report(optim.sys.output,'start krotov optimisation using quasi-newton, Symmetric Rank 1 Hessian update',1);
            case 'krotov-bfgs', optim_report(optim.sys.output,'start krotov optimisation using quasi-newton, BFGS method',1);
        end
        optim_report(optim.sys.output,['Cost Function = @', func2str(optim.sys.costfun.handle)],1);
        optim_report(optim.sys.output,'==============================================================================================',1);
        obj_fun_str=[blanks(17) optim.sys.obj_fun_str]; obj_fun_str=obj_fun_str(end-17:end);
        if ismember(optim.reg.method,{'RFO','TRM'}) || optim.reg.track_eig
            optim_report(optim.sys.output,['Iter  Time #f(x) #g(x)  min(eig) #cond  cond(H) ' obj_fun_str '  1st-optim'],1);
        elseif ismember(optim.reg.method,{'CHOL','reset'})
            optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #reg ' obj_fun_str '  1st-optim'],1);
        else optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) ' obj_fun_str '  1st-optim'],1);
        end
        optim_report(optim.sys.output,'----------------------------------------------------------------------------------------------',1);
    case {'iter'}
        t=data.itertime;
        if t<1e-3, Tstr=sprintf('%5.0e',t); elseif  t<10, Tstr=sprintf('%5.3f',t); elseif  t<1e2, Tstr=sprintf('%5.2f', t);
        elseif  t<1e3, Tstr=sprintf('%5.1f',t); elseif  t<=99999, Tstr=sprintf('%5.0f',t); else Tstr=sprintf('%5.0e',t);
        end
        % Report optimisation iterasation to user
        if (ismember(optim.reg.method,{'RFO','TRM'}) && (isfield(data,'n_cond') && data.n_cond>0)) || optim.reg.track_eig
            optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f  %8.3g %4.0f %9.3g  %17.9g %9.3g'],...
                data.iter,data.n_fx,data.n_grad,data.hess_mineig,data.n_cond,data.hess_cond,optim.sys.obj_fun_disp(fx),data.firstorderopt),1);
        elseif ismember(optim.reg.method,{'RFO','TRM'})
            optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f  %8.1s %4.0f %9.1s  %17.9g %9.3g'],...
                data.iter,data.n_fx,data.n_grad,'-',0,'-',optim.sys.obj_fun_disp(fx),data.firstorderopt),1);
        elseif ismember(optim.reg.method,{'CHOL','reset'})
            optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %17.9g %9.3g'],...
                data.iter,data.n_fx,data.n_grad,data.n_cond,optim.sys.obj_fun_disp(fx),data.firstorderopt),1);
        else optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f  %17.9g %9.3g'],...
                data.iter,data.n_fx,data.n_grad,optim.sys.obj_fun_disp(fx),data.firstorderopt),1);
        end
        
    case {'term'}
        optim_report(optim.sys.output,'==============================================================================================',1);
        % Make exist output structure
        switch(optim.algo.method)
            case 'krotov', data.algorithm='Krotov method';
            case 'krotov-sr1', data.algorithm='Krotov method, Symmetric Rank 1 Hessian update (SR1)';
            case 'krotov-bfgs', data.algorithm='Krotov method, Broyden-Fletcher-Goldfarb-Shanno Hessian update (BFGS)';
        end
        
        % report optimisation results to user
        data.exit_message=exitmessage(exitflag);
        data.minimum = optim.sys.obj_fun_disp(fx);
        data.timeTotal=toc(data.timeTotal);
        data.timeIntern=data.timeTotal-data.time_fx-data.time_gfx;
        optim_report(optim.sys.output,'Optimisation Results',1)
        optim_report(optim.sys.output,['    Algorithm Used     : ' data.algorithm],1);
        optim_report(optim.sys.output,['    Exit message       : ' data.exit_message],1);
        optim_report(optim.sys.output,['    Iterations         : ' int2str(data.iter)],1);
        optim_report(optim.sys.output,['    Function Count     : ' int2str(data.n_fx)],1);
        optim_report(optim.sys.output,['    Gradient Count     : ' int2str(data.n_grad)],1);
        optim_report(optim.sys.output,['    Minimum found      : ' num2str(data.minimum)],1);
        optim_report(optim.sys.output,['    Optimiser Time     : ' num2str(data.timeIntern) ' seconds'],1);
        optim_report(optim.sys.output,['    Function calc Time : ' num2str(data.time_fx) ' seconds'],1);
        optim_report(optim.sys.output,['    Gradient calc Time : ' num2str(data.time_gfx) ' seconds'],1);
        optim_report(optim.sys.output,['    Total Time         : ' num2str(data.timeTotal) ' seconds'],1);
        optim_report(optim.sys.output,'==============================================================================================',1);
        if ~isinteger(optim.sys.output) && optim.sys.output>=3, optim_report(optim.sys.output,[],[],true); end
        pause(1.5);
end
end

function message=exitmessage(exitflag)
switch(exitflag) % list of exiflag messages
    case  1, message='Objective function, f(x), is not strictly decreasing';
    case  2, message='Change in x was smaller than the specified tolerance, tol_x.';
    case  3, message='Boundary fminimum reached.';
    case  4, message='Specified terminal functional value reached';
    case  0, message='Number of iterations exceeded the specified maximum.';
    case -1, message='Algorithm was terminated by the output function.';
    case -2, message='Line search cannot find an acceptable point along the current search';
    case -3, message='Hessian matrix is not positive definite';
    otherwise, message='Undefined exit code';
end
end

% There is a sacred horror about everything grand. It is easy to admire
% mediocrity and hills; but whatever is too lofty, a genius as well as a
% mountain, an assembly as well as a masterpiece, seen too near, is
% appalling... People have a strange feeling of aversion to anything grand.
% They see abysses, they do not see sublimity; they see the monster, they
% do not see the prodigy.
%
% Victor Hugo - Ninety-three

