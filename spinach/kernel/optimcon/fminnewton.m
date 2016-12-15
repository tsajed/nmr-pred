% Finds a local minimum of a function of several variables using
% newton-based algorithms.
%
% Based on fminlbfgs.m code from D. Kroon, University of Twente (Nov 2010).
%
% <http://spindynamics.org/wiki/index.php?title=Fminnewton.m>

function [x,fx,grad,hess,data]=fminnewton(cost_function,x_0,optim,hess_init,cost_fun_vars)

% Initialize the structure of the initial guess
x = x_0(:); data.x_shape=size(x_0); data.n_vars = numel(x_0);

% Always Bootstrap and set optimisation tolerances
% - assume user isn't optimisation expert
optim.cost_function=cost_function;
optim=optim_tols(optim,data.n_vars);

% initialise counters and timers
data.iter=0; data.n_fx=0; data.n_grad=0; data.n_hess=0; data.n_cond=0; data.firstorderopt=[];
data.timeTotal=tic; data.time_fx=0; data.time_gfx=0; data.time_hfx=0; data.iter_timer=tic;

% add cost function variables to the data structure (empty if not provided)
if ~exist('cost_fun_vars','var'), cost_fun_vars=[]; end
data.cost_fun_vars=cost_fun_vars;

% penalty terms to consider (default=0), affects output display only.
data.npen=optim.algo.npen;

% initialise lbfgs memory counter
if ismember(optim.algo.method,{'lbfgs'}),
    store.lbfgs_store=0;
    store.n_store=optim.algo.n_store; optim.algo=rmfield(optim.algo,'n_store');
    store.deltaX=optim.algo.deltaX; optim.algo=rmfield(optim.algo,'deltaX');
    store.deltaG=optim.algo.deltaG; optim.algo=rmfield(optim.algo,'deltaG');
end

% print optimisation header
iterative_reporting([],[],'head',[],optim,data);

% set initial Hessian approximation if given, for bfgs or sr1
if ismember(optim.algo.method,{'bfgs','sr1'}) && exist('hess_init','var')
    hess=hess_init;
else
    hess=eye(data.n_vars);
end

% initialise exitflag, functional and its derivatives
exitflag=[]; fx=[]; grad=[]; alpha=0; dir=[];

% send to output function if defined
if(output_function(optim.sys.outputfun,'init')), exitflag=-1; end

% Start optimisation
while(true)
    
    % Calculate the initial error and gradient and hessian
    if ismember(optim.algo.method,{'newton'}) && isempty(exitflag)
        [data,fx,grad,hess]=objeval(x(:),cost_function, data);
    elseif isempty(grad) && isempty(exitflag)
        [data,fx,grad]=objeval(x(:),cost_function, data);
    end
    
    % check if iterations exceeds maximum
    if(optim.tols.max_iter <= data.iter), exitflag=0; end
    
    %search direction (p_current)
    if(isempty(exitflag)), [dir,grad,hess,data]=search_dir(grad,hess,optim,data,dir); end
    
    if(isempty(exitflag)), data.firstorderopt=norm(grad,Inf); end
    
    % Show the current iteration information
    data.itertime=toc(data.iter_timer); data.iter_timer=tic;
    
    % report diagnostics to user
    iterative_reporting(fx,alpha,'iter',exitflag,optim,data);
    
    % Call output function if specified
    if(output_function(optim.sys.outputfun,'iter')), exitflag=-1; end
    
    % Check if exitflag is set
    if(~isempty(exitflag)), break, end
    
    % Keep the variables for next iteration
    store.fx=fx; store.x=x(:); store.grad=grad(:);
    
    % make a linesearch in the search direction
    [alpha,fx,grad,exitflag,data] = linesearch(cost_function,...
        dir, x, fx, grad, data, optim);
    if isempty(fx), fx=store.fx; grad=store.grad; x=store.x; end
    
    % update current x with search direction and step length
    x=x+alpha*dir;
    
    % hessian update if quasi-newton
    if ismember(optim.algo.method,{'lbfgs','bfgs','sr1'}) && (isempty(exitflag))
        % required when the linesearch has used the newton step
        if isempty(grad)
            [data,fx,grad]=objeval(x(:),cost_function, data);
        end
        [hess,dir,store]=hessian_update(x,grad,hess,store,optim.algo);
    end
    
    % Update number of iterations
    if (isempty(exitflag)), data.iter=data.iter+1; end
    
    % Gradient norm too small, change in x too small, terminal functional value reached.
    if(optim.algo.max_min*optim.tols.fx  < optim.algo.max_min*fx), exitflag=4; end
    if(optim.tols.x   > max(abs(store.x-x))), exitflag=2; end
    if(optim.tols.gfx > data.firstorderopt), exitflag=1; end
    
end

% Call output function if specified
if(output_function(optim.sys.outputfun,'done')), exitflag=-1; end

% ensure outputs exist, if not then re-evaluate
if ismember(optim.algo.method,{'newton'}) && isempty(grad)
    [data,fx,grad,hess]=objeval(x(:),cost_function, data);
elseif isempty(grad)
    [data,fx,grad]=objeval(x(:),cost_function, data);
end

% Set outputs
x=reshape(x,data.x_shape); grad=reshape(grad,data.x_shape);
data=iterative_reporting(fx,alpha,'term',exitflag,optim,data);

    function exit=output_function(outputfun,where)
        exit=false; % initialise
        if(~isempty(outputfun)) % data structure to output function
            % <add any variable to export here>
            % if exit is returned as true, optimisation will finish
            exit=feval(outputfun,reshape(x,data.x_shape),data,where);
        end
    end
end

function [dir,grad,hess,data]=search_dir(grad,hess,optim,data,dir)

% find the search direction
switch optim.algo.method
    case {'grad_descent','grad_ascent'}
        
        % set the search direction as the gradient direction
        dir=grad;
        
    case 'lbfgs'
        
        % calculate the direction for the first iterate only
        if data.iter==0, dir=optim.algo.max_min.*grad;  end
        
    case {'sr1','bfgs','newton'}
        
        % prepare the Hessian to be definite and well-conditioned
        [grad,hess,data]=hessprep(grad,hess,optim,data);
        
        % force real and symmetric hessian
        if ~isreal(hess), hess=real(hess); grad=real(grad); end
        if ~issymmetric(hess), hess=(hess+hess')./2; end
        
        % calculate the search direction from the Hessian and the gradient
        if xor((strcmp(optim.reg.method,'CHOL') && data.n_cond>0), optim.algo.inv_method)
            dir = -hess*grad;
        else dir = -hess\grad;
        end
end

end

function [hess,dir,store] = hessian_update(x,grad,hess,store,algo)

switch algo.method
    
    case 'lbfgs' % limited-memory-BFGS
        
        % Update a list with the history of deltaX and deltaG
        store.deltaX=[x-store.x         store.deltaX(:,1:store.n_store-1)];
        store.deltaG=[grad-store.grad	store.deltaG(:,1:store.n_store-1)];
        
        % size of store of deltaX and deltaG
        store.lbfgs_store=min(store.lbfgs_store+1, store.n_store);
        
        % Set new Direction with L-BFGS with scaling as described by Nocedal
        dir=lbfgs(store.deltaX, store.deltaG, grad, store.lbfgs_store);
        
    case {'bfgs','sr1'}
        
        % calculate BFGS Hessian update
        hess=quasinewton(hess, x, store.x, -algo.max_min*grad, -algo.max_min*store.grad, algo.inv_method, algo.method);
        dir=[];     % no search direction yet
        
    otherwise
        error('unsupported quasi-Newton hessian update method')
        
end

end

function data=iterative_reporting(fx,alpha,report_case,exitflag,optim,data)

switch report_case
    case {'head'}
        % Report start of optimisation to user
        switch(optim.algo.method)
            case 'sr1', optim_report(optim.sys.output,['start ' optim.algo.extremum(1:5) 'isation using quasi-newton, Symmetric Rank 1 Hessian update'],1);
            case 'bfgs', optim_report(optim.sys.output,['start ' optim.algo.extremum(1:5) 'isation using quasi-newton, BFGS method'],1);
            case 'lbfgs',  optim_report(optim.sys.output,['start ' optim.algo.extremum(1:5) 'isation using quasi-newton, limited memory BFGS method'],1);
            case 'grad_descent', optim_report(optim.sys.output,['start ' optim.algo.extremum(1:5) 'isation using Steepest Gradient Descent method'],1);
            case 'grad_ascent', optim_report(optim.sys.output,['start ' optim.algo.extremum(1:5) 'isation using Steepest Gradient Ascent method'],1);
            case 'newton', optim_report(optim.sys.output,['start ' optim.algo.extremum(1:5) 'isation using Newton-Raphson method'],1);
        end
        optim_report(optim.sys.output,['Cost Function = @', func2str(optim.sys.costfun.handle)],1);
        obj_fun_str=[blanks(17) optim.sys.obj_fun_str]; obj_fun_str=obj_fun_str(end-17:end);
        if ~optim.algo.pen
            optim_report(optim.sys.output,'==============================================================================================',1);
            if ismember(optim.reg.method,{'RFO','TRM'}) || optim.reg.track_eig
                optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #H(x)  min|eig| #cond  cond(H) ' obj_fun_str '   Step-size 1st-optim'],1);
            elseif ismember(optim.reg.method,{'CHOL','reset'})
                optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #H(x) #reg ' obj_fun_str '   Step-size 1st-optim'],1);
            else optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #H(x) ' obj_fun_str '   Step-size 1st-optim'],1);
            end
            optim_report(optim.sys.output,'----------------------------------------------------------------------------------------------',1);
        else
            optim_report(optim.sys.output,'==========================================================================================================',1);
            if ismember(optim.reg.method,{'RFO','TRM'}) || optim.reg.track_eig
                optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #H(x)  min|eig| #cond  cond(H) ' obj_fun_str '    sum(pen)   Step-size 1st-optim'],1);
            elseif ismember(optim.reg.method,{'CHOL','reset'})
                optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #H(x) #reg ' obj_fun_str '    sum(pen)   Step-size 1st-optim'],1);
            else optim_report(optim.sys.output,['Iter  Time #f(x) #g(x) #H(x) ' obj_fun_str '    sum(pen)   Step-size 1st-optim'],1);
            end
            optim_report(optim.sys.output,'----------------------------------------------------------------------------------------------------------',1);
        end
    case {'iter'}
        % Report optimisation iterate to user
        t=data.itertime;
        if t<1e-3, Tstr=sprintf('%5.0e',t); elseif  t<10, Tstr=sprintf('%5.3f',t); elseif  t<1e2, Tstr=sprintf('%5.2f', t);
        elseif  t<1e3, Tstr=sprintf('%5.1f',t); elseif  t<=99999, Tstr=sprintf('%5.0f',t); else Tstr=sprintf('%5.0e',t);
        end
        if ~optim.algo.pen
            if (ismember(optim.reg.method,{'RFO','TRM'}) && (isfield(data,'n_cond') && data.n_cond>0)) || optim.reg.track_eig
                optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %8.3g %4.0f %9.3g  %17.9g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,data.hess_mineig,data.n_cond,data.hess_cond,optim.sys.obj_fun_disp(fx),alpha,data.firstorderopt),1);
            elseif ismember(optim.reg.method,{'RFO','TRM'})
                optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %8.1s %4.0f %9.1s  %17.9g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,'-',0,'-',optim.sys.obj_fun_disp(fx),alpha,data.firstorderopt),1);
            elseif ismember(optim.reg.method,{'CHOL','reset'})
                optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f %5.0f  %17.9g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,data.n_cond,optim.sys.obj_fun_disp(fx),alpha,data.firstorderopt),1);
            else optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %17.9g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,optim.sys.obj_fun_disp(fx),alpha,data.firstorderopt),1);
            end
        else
            if (ismember(optim.reg.method,{'RFO','TRM'}) && (isfield(data,'n_cond') && data.n_cond>0)) || optim.reg.track_eig
                optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %8.3g %4.0f %9.3g  %17.9g  %10.4g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,data.hess_mineig,data.n_cond,data.hess_cond,...
                    optim.sys.obj_fun_disp(data.fx_sep_pen(1)),optim.sys.obj_fun_disp(sum(data.fx_sep_pen(2:end))),alpha,data.firstorderopt),1);
            elseif ismember(optim.reg.method,{'RFO','TRM'})
                optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %8.1s %4.0f %9.1s  %17.9g  %10.4g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,'-',0,'-',...
                    optim.sys.obj_fun_disp(data.fx_sep_pen(1)),optim.sys.obj_fun_disp(sum(data.fx_sep_pen(2:end))),alpha,data.firstorderopt),1);
            elseif ismember(optim.reg.method,{'CHOL','reset'})
                optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f %5.0f  %17.9g  %10.4g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,data.n_cond,...
                    optim.sys.obj_fun_disp(data.fx_sep_pen(1)),optim.sys.obj_fun_disp(sum(data.fx_sep_pen(2:end))),alpha,data.firstorderopt),1);
            else optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f %5.0f %5.0f  %17.9g  %10.4g  %10.5g %9.3g'],...
                    data.iter,data.n_fx,data.n_grad,data.n_hess,...
                    optim.sys.obj_fun_disp(data.fx_sep_pen(1)),optim.sys.obj_fun_disp(sum(data.fx_sep_pen(2:end))),alpha,data.firstorderopt),1);
            end
        end
        
    case {'term'}
        % Make exit output structure
        optim_report(optim.sys.output,'==============================================================================================',1);
        switch(optim.algo.method)
            case 'sr1', data.algorithm='Symmetric Rank 1 Hessian update (SR1)';
            case 'bfgs', data.algorithm='Broyden-Fletcher-Goldfarb-Shanno Hessian update (BFGS)';
            case 'lbfgs',  data.algorithm='limited memory Broyden-Fletcher-Goldfarb-Shanno update (L-BFGS)';
            case 'grad_descent', data.algorithm='Steepest Gradient Descent';
            case 'grad_ascent', data.algorithm='Steepest Gradient Ascent';
            case 'newton', data.algorithm='Newton-Raphson method';
        end
        
        % report optimisation results to user
        data.exit_message=exitmessage(exitflag);
        data.extremum = optim.sys.obj_fun_disp(fx);
        data.timeTotal=toc(data.timeTotal);
        data.timeIntern=data.timeTotal-data.time_fx-data.time_gfx-data.time_hfx;
        optim_report(optim.sys.output,'Optimisation Results',1)
        optim_report(optim.sys.output,['    Algorithm Used     : ' data.algorithm],1);
        optim_report(optim.sys.output,['    Exit message       : ' data.exit_message],1);
        optim_report(optim.sys.output,['    Iterations         : ' int2str(data.iter)],1);
        optim_report(optim.sys.output,['    Function Count     : ' int2str(data.n_fx)],1);
        optim_report(optim.sys.output,['    Gradient Count     : ' int2str(data.n_grad)],1);
        optim_report(optim.sys.output,['    Hessian Count      : ' int2str(data.n_hess)],1);
        optim_report(optim.sys.output,['    M' optim.algo.extremum(2:end) ' found      : ' num2str(data.extremum)],1);
        optim_report(optim.sys.output,['    norm(gradient)     : ' num2str(data.firstorderopt)],1);
        optim_report(optim.sys.output,['    Optimiser Time     : ' num2str(data.timeIntern) ' seconds'],1);
        optim_report(optim.sys.output,['    Function calc Time : ' num2str(data.time_fx) ' seconds'],1);
        optim_report(optim.sys.output,['    Gradient calc Time : ' num2str(data.time_gfx) ' seconds'],1);
        optim_report(optim.sys.output,['    Hessian calc Time  : ' num2str(data.time_hfx) ' seconds'],1);
        optim_report(optim.sys.output,['    Total Time         : ' num2str(data.timeTotal) ' seconds'],1);
        optim_report(optim.sys.output,'==============================================================================================',1);
        if isnumeric(optim.sys.output) && optim.sys.output>=3, optim_report(optim.sys.output,[],[],true); end
end
end

function message=exitmessage(exitflag)
% list of exiflag messages
switch(exitflag)
    case  1, message='Change in gradient norm less than specified tolerance, tol_gfx.';
    case  2, message='Change in x was smaller than the specified tolerance, tol_x.';
    case  3, message='Boundary reached.';
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

