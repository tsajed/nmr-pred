% Nelder-Mead and Multi-directional search simplex algorithms for
% minimisation of a gradient-free objective function
%
% <http://spindynamics.org/wiki/index.php?title=Fminsimplex.m>

function [x, f_min, data] = fminsimplex(cost_function, x, optim, cost_fun_vars)

x0 = x(:); data.x_shape=size(x); % Work with column vector internally.
data.n_vars = numel(x0);

% Always Bootstrap and set optimisation tolerances
% - assume user isn't optimisation expert
optim.cost_function=cost_function;
optim=optim_tols(optim,data.n_vars);

% add cost function variables to the data structure (empty if not provided)
if ~exist('cost_fun_vars','var'), cost_fun_vars=[]; end
data.cost_fun_vars=cost_fun_vars;

% penalty terms to consider (default=0), affects output display only.
data.npen=optim.algo.npen;

X = [zeros(data.n_vars,1) eye(data.n_vars)];

X(:,1) = x0; data.n_fx=0;
data.iter = 0; data.iter_inner = 0;
data.parallel_cost_fcn=optim.algo.parallel_cost_fcn;
data.timeTotal=tic; data.time_fx=0; data.iter_timer=tic;

% Set up initial simplex.
scale = max(norm(x0,inf),1);
switch optim.algo.init_simplex
    case {'equilateral'}
        % Regular simplex - all edges have same length.
        % Generated from construction given in reference [18, pp. 80-81] of [1].
        alpha = scale / (data.n_vars*sqrt(2)) * [ sqrt(data.n_vars+1)-1+data.n_vars  sqrt(data.n_vars+1)-1 ];
        X(:,2:data.n_vars+1) = (x0 + alpha(2)*ones(data.n_vars,1)) * ones(1,data.n_vars);
        for j=2:data.n_vars+1
            X(j-1,j) = x0(j-1) + alpha(1);
        end
    case {'right-angled'}
        % Right-angled simplex based on co-ordinate axes.
        alpha = scale*ones(data.n_vars+1,1);
        for j=2:data.n_vars+1
            X(:,j) = x0 + alpha(j)*X(:,j);
        end
end
% find simplex functional values (cost function should be parallel)
[data,f_x]=objeval(X,cost_function, data);
fmin_old = f_x(1);

% print optimisation header
iterative_reporting([],[],'head',[],optim,data);
exitflag=[]; data.how=' initial'; data.size = 0;
data.n_fx = data.n_vars+1;
[f_x,j] = sort(f_x); X = X(:,j);
f_min = f_x(1); x_min = X(:,1);
data.size_simplex = norm(X(:,2:data.n_vars+1)-x_min*ones(1,data.n_vars),1) / max(1, norm(x_min,1));

% Start optimisation
while(true)
    
    % define the minimum
    f_min = f_x(1);
    
    % Add one to the current iteration
    data.itertime=toc(data.iter_timer); data.iter_timer=tic; data.iter = data.iter+1;
    
    % check if iterations exceeds maximum
    if(optim.tols.max_iter < data.iter), exitflag=-1; break, end
    
    if f_min < fmin_old
        % report diagnostics to user
        data.size_simplex = norm(X(:,2:data.n_vars+1)-x_min*ones(1,data.n_vars),1) / max(1, norm(x_min,1));
        iterative_reporting(f_min,fmin_old,'iter',exitflag,optim,data);
        
        % Call output function if specified
        if(output_function(optim.sys.outputfun,'iter')), exitflag=-1; break, end
    end
    x_min = X(:,1); fmin_old = f_min;
    
    % Stopping Test 1 - f reached target value?
    if f_min <= optim.tols.termination ,exitflag = 0; break, end
    
    switch optim.algo.method
        case 'md_simplex'
            data.iter_inner = 0;
            while(true)   %%% Inner repeat loop.
                data.iter_inner = data.iter_inner+1;
                
                % Stopping Test 2 - too many f-evals?
                if data.n_fx >= optim.tols.max_n_fx, exitflag = 1; break, end
                
                % Stopping Test 3 - converged?   test (4.3) in [1].
                data.size_simplex = norm(X(:,2:data.n_vars+1)-x_min*ones(1,data.n_vars),1) / max(1, norm(x_min,1));
                
                if data.size_simplex <= optim.tols.simplex_min, exitflag = 2; break, end
                
                X_r=zeros(length(x_min),data.n_vars);
                for j=1:data.n_vars      % ---Rotation (reflection) step.
                    X_r(:,j) = x_min - (X(:,j+1)-x_min);
                end
                [data,f_r]=objeval(X_r,cost_function, data);
                
                replaced = ( min(f_r) < f_min);
                
                if replaced
                    X_e=zeros(length(x_min),data.n_vars);
                    for j=1:data.n_vars   % ---Expansion step.
                        X_e(:,j) = x_min - optim.tols.expansion*(X_r(:,j)-x_min);
                    end
                    [data,f_e]=objeval(X_e,cost_function, data);
                    
                    % Accept expansion or rotation?
                    if min(f_e) < min(f_x)
                        X(:,2:data.n_vars+1) = X_e;
                        f_x(2:data.n_vars+1) = f_e;  % Accept rotation.
                        data.size = data.size + 1;  % Accept expansion (f and V already set).
                        data.how='  expand';
                    else
                        X(:,2:data.n_vars+1) = X_r;
                        f_x(2:data.n_vars+1) = f_r;  % Accept rotation.
                        data.how=' reflect';
                    end
                else
                    X_c=zeros(length(x_min),data.n_vars);
                    for j=1:data.n_vars   % ---Contraction step.
                        X_c(:,j) = x_min + optim.tols.contraction*(X(:,j+1)-x_min);
                    end
                    [data,f_c]=objeval(X_c,cost_function, data);
                    
                    replaced = (min(f_c) < f_min);
                    data.how='contract';
                    % Accept contraction (f and V already set).
                    
                    X(:,2:data.n_vars+1) = X_c;
                    f_x(2:data.n_vars+1) = f_c;  % Accept rotation.
                    data.size = data.size - 1;
                end
                
                [f_x,j]=sort(f_x); X = X(:,j);
                
                if replaced, break, end
                if optim.sys.output && rem(data.iter_inner,10) == 0
                    optim_report(optim.sys.output,['        ...inner iterations = ' int2str(data.iter_inner)],1);
                end
            end %%% Of inner repeat loop.
            
            if ~isempty(exitflag), break, end
        case 'nm_simplex'
            % Stopping Test 2 - too many f-evals?
            if data.n_fx >= optim.tols.max_n_fx, exitflag = 1; break, end
            
            % Stopping Test 3 - converged?   test (4.3) in [1].
            x_min = X(:,1);
            data.size_simplex = norm(X(:,2:data.n_vars+1)-x_min*ones(1,data.n_vars),1) / max(1, norm(x_min,1));
            
            if data.size_simplex <= optim.tols.simplex_min, exitflag = 2; break, end
            
            %  One step of the Nelder-Mead simplex algorithm
            vbar = sum(X(:,1:data.n_vars),2)./data.n_vars;  % Mean value
            X_r = (1 + optim.tols.reflection)*vbar - optim.tols.reflection*X(:,data.n_vars+1); x = X_r(:);
            
            [data,f_r]=objeval(x,cost_function, data);
            X_k = X_r;  f_k = f_r; data.how = ' reflect';
            if f_r < f_x(1)
                if f_r < f_x(1)
                    X_e = optim.tols.expansion*X_r + (1-optim.tols.expansion)*vbar; x = X_e(:);
                    [data,f_e]=objeval(x,cost_function, data);
                    if f_e < f_x(1)
                        X_k = X_e; f_k = f_e;
                        data.how = '  expand';
                    end
                end
            else
                X_t = X(:,1); f_t = f_x(1);
                if f_r < f_t
                    X_t = X_r;  f_t = f_r;
                end
                X_c = optim.tols.contraction*X_t + (1-optim.tols.contraction)*vbar; x = X_c(:);
                [data,f_c]=objeval(x,cost_function, data);
                if f_c < f_x(data.n_vars)
                    X_k = X_c; f_k = f_c;
                    data.how = 'contract';
                else
                    for j = 2:data.n_vars
                        X(:,j) = (X(:,1) + X(:,j))/2;
                    end
                    [data,f_x]=objeval(X,cost_function, data);
                    
                    X_k = (X(:,1) + X(:,data.n_vars+1))/2; x(:) = X_k;
                    [data,f_k]=objeval(x,cost_function, data);
                    data.how = '  shrink';
                end
            end
            X(:,data.n_vars+1) = X_k;
            f_x(data.n_vars+1) = f_k;
            [f_x,j] = sort(f_x); X = X(:,j);
    end
end   %%%%%% Of outer loop.

% Call output function if specified
if(output_function(optim.sys.outputfun,'done')), exitflag=-2; end

% Finished.
switch optim.algo.method
    case 'md_simplex'
        x(:) = x_min;
        x=reshape(x,data.x_shape);
    case 'nm_simplex'
        x(:) = X(:,1);
        x=reshape(x,data.x_shape);
end
data=iterative_reporting(f_min,fmin_old,'term',exitflag,optim,data);

    function exit=output_function(outputfun,where)
        exit=false; % initialise
        if(~isempty(outputfun)) % data structure to output function
            % <add any variable to export here>
            % if exit is returned as true, optimisation will finish
            if isempty(opt_arg)
                exit=feval(outputfun,reshape(x,data.x_shape),data,where);
            else
                exit=feval(outputfun,reshape(x,data.x_shape),data,where,opt_arg);
            end
        end
    end
end

function data=iterative_reporting(fx_min,fmin_old,report_case,exitflag,optim,data)

switch report_case
    case {'head'}
        % Report start of optimisation to user
        switch(optim.algo.method)
            case 'nm_simplex', optim_report(optim.sys.output,'start optimisation using Nelder-Mead, simplex method',1);
            case 'md_simplex', optim_report(optim.sys.output,'start optimisation using Multi-directional, simplex method',1);
        end
        optim_report(optim.sys.output,['Cost Function = @', func2str(optim.sys.costfun.handle)],1);
        optim_report(optim.sys.output,'==============================================================================================',1);
        obj_fun_str=[blanks(17) optim.sys.obj_fun_str]; obj_fun_str=obj_fun_str(end-17:end);
        if ismember(optim.algo.method,{'nm_simplex'})
            optim_report(optim.sys.output,['Iter  Time #f(x)      how  ' obj_fun_str ' splx-size   Change'],1);
        elseif ismember(optim.algo.method,{'md_simplex'})
            optim_report(optim.sys.output,['Iter  Time #f(x)      how InIter ' obj_fun_str ' splx-size   Change'],1);
        end
        optim_report(optim.sys.output,'----------------------------------------------------------------------------------------------',1);
    case {'iter'}
        t=data.itertime;
        if t<1e-3, Tstr=sprintf('%5.0e',t); elseif  t<10, Tstr=sprintf('%5.3f',t); elseif  t<1e2, Tstr=sprintf('%5.2f', t);
        elseif  t<1e3, Tstr=sprintf('%5.1f',t); elseif  t<=99999, Tstr=sprintf('%5.0f',t); else Tstr=sprintf('%5.0e',t);
        end
        % Report optimisation iterasation to user
        if ismember(optim.algo.method,{'nm_simplex'})
            optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f ' data.how ' %17.9g   %10.5g  (%2.1f%%)'],...
                data.iter,data.n_fx,optim.sys.obj_fun_disp(fx_min),data.size_simplex,100*(fx_min-fmin_old)/(abs(fmin_old)+eps)),1);
        elseif ismember(optim.algo.method,{'md_simplex'})
            optim_report(optim.sys.output,sprintf(['%4.0f ',Tstr,' %5.0f ' data.how ' %5.0f %17.9g   %10.5g  (%2.1f%%)'],...
                data.iter,data.n_fx,data.iter_inner,optim.sys.obj_fun_disp(fx_min),data.size_simplex,100*(fx_min-fmin_old)/(abs(fmin_old)+eps)),1);
        end
    case {'term'}
        optim_report(optim.sys.output,'==============================================================================================',1);
        % Make exist output structure
        switch(optim.algo.method)
            case 'nm_simplex', data.algorithm='Nelder-Mead Simplex method';
            case 'md_simplex', data.algorithm='Multi-directional Simplex method';
        end
        
        % report optimisation results to user
        data.exit_message=exitmessage(exitflag);
        data.minimum = optim.sys.obj_fun_disp(fx_min);
        data.timeTotal=toc(data.timeTotal);
        data.timeIntern=data.timeTotal-data.time_fx;
        optim_report(optim.sys.output,'Optimisation Results',1)
        optim_report(optim.sys.output,['    Algorithm Used     : ' data.algorithm],1);
        optim_report(optim.sys.output,['    Exit message       : ' data.exit_message],1);
        optim_report(optim.sys.output,['    Iterations         : ' int2str(data.iter)],1);
        optim_report(optim.sys.output,['    Function Count     : ' int2str(data.n_fx)],1);
        optim_report(optim.sys.output,['    Minimum found      : ' num2str(data.minimum)],1);
        optim_report(optim.sys.output,['    Optimiser Time     : ' num2str(data.timeIntern) ' seconds'],1);
        optim_report(optim.sys.output,['    Function calc Time : ' num2str(data.time_fx) ' seconds'],1);
        optim_report(optim.sys.output,['    Total Time         : ' num2str(data.timeTotal) ' seconds'],1);
        optim_report(optim.sys.output,'==============================================================================================',1);
        if ~isinteger(optim.sys.output) && optim.sys.output>=3, optim_report(optim.sys.output,[],[],true); end
        pause(1.5);
end
end

function message=exitmessage(exitflag)
switch(exitflag) % list of exiflag messages
    case  2, message='Simplex size less than defined tolerance <simplex_min>';
    case  1, message='Max no. of function evaluations exceeded.';
    case  0, message='Exceeded target.';
    case -1, message='Max no. of iterations exceeded.';
    case -2, message='Algorithm was terminated by the output function.';
    otherwise, message='Undefined exit code';
end
end

% Seek freedom and become captive of your desires. Seek discipline and find
% your liberty.
%
% Frank Herbert
