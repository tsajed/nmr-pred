% Tolerances and options for optimisation and control. Sets the various
% optimisation methods, corresponding tolerances, and options used
% throughout numerical optimisation.
%
% Modifications to this function are discouraged -- the accuracy settings
% should be modified by setting the sub-fields of the optim structure.
%
% <http://spindynamics.org/wiki/index.php?title=Optim_tols.m>

function optim_param=optim_tols(optim,n_vars)

% list of OPTIMISATION algorithms
% >>>>-add to list as appropriate-<<<<
optim_algos={'lbfgs','bfgs','newton','sr1','grad_ascent','grad_descent','krotov','krotov-bfgs','krotov-sr1'};

% list of LINE-SEARCH methods and rules
% >>>>-add to lists as appropriate-<<<<
search_algos={'bracket-section','backtracking','newton-step'};
search_rules={'Wolfe','Wolfe-strong','Wolfe-weak','Goldstein','Armijo'};

% list of REGULARISATION and conditioning proceedures
% >>>>-add to lists as appropriate-<<<<
reg_method={'RFO','TRM','CHOL'};
reg_cond={'iterative','scaled','none'};

% max number of variables before switch to lbfgs (excluding grad_descent and grad_ascent)
lbfgs_switch=5000;

% ========================== STRUCT = OPTIM.SYS ==========================
% store the cost function handle; throw error is not supplied - should
% always be supplied if used with fminnewton: inherited from fminnewton
% inputs
if isempty(optim)
    error('Empty optimisation options');
elseif ~isfield(optim,'cost_function') || ~isa(optim.cost_function,'function_handle'),
    error('Optimisation functional must be supplied as a function handle e.g. @cost_function')
else
    optim_param.sys.costfun=functions(optim.cost_function);
    optim_param.sys.costfun.handle=optim.cost_function;
    optim=rmfield(optim,'cost_function'); % parsed -> remove
    % remove workspace field - can be commented if appropriate
    if isfield(optim_param.sys.costfun,'workspace')
        optim_param.sys.costfun = rmfield(optim_param.sys.costfun,'workspace');
    end
end

% Decide output destination
if isfield(optim,'output') && strcmp(optim.output,'hush')
    % Hush the output
    optim_param.sys.output='hush';
    optim=rmfield(optim,'output'); % parsed -> remove
elseif isfield(optim,'output')&&(~strcmp(optim.output,'console'));
    % Close all open files
    
    % Print to a user-specified file
    optim_param.sys.output=fopen(optim.output,'a');
    optim=rmfield(optim,'output'); % parsed -> remove
elseif isfield(optim,'output')&&(strcmp(optim.output,'console'));
    optim_param.sys.output=1;
    optim=rmfield(optim,'output'); % parsed -> remove
else
    % Print to the console (default -> force)
    optim_param.sys.output=1;
end

% Optimisation output function, evaluated at each iteration
if isfield(optim,'outputfun') && isa(optim.outputfun,'function_handle')
    optim_param.sys.outputfun=optim.outputfun;
    optim=rmfield(optim,'outputfun'); % parsed -> remove
else
    % (default -> force)
    optim_param.sys.outputfun='';
end
% get function call info for cost_function
file_parse_info=getcallinfo(optim_param.sys.costfun.file); %optim_param.sys.costfun.info
% ========================================================================

if ismember(optim_param.sys.costfun.type,{'scopedfunction'})
    optim_param.sys.costfun.type='local'; %scopedfunction=local function
end
if ismember(optim_param.sys.costfun.type,{'nested','local','simple'})
    [~,optim_param.sys.costfun.name]=fileparts(optim_param.sys.costfun.function);
    for n=1:numel(file_parse_info)
        if strcmp(file_parse_info(n).name,optim_param.sys.costfun.name)
            optim_param.sys.costfun.fcn_calls=unique(file_parse_info(n).calls.fcnCalls.names)';
        end
    end
elseif strcmp(optim_param.sys.costfun.type,'anonymous')
    % throw error for anonymous function
    error('cost function handle should not be anonymous - please use nested, local, or simple functions')
end

fcn_list=dbstack;
if numel(fcn_list)==1, optimise_fcn={}; end
if numel(fcn_list)>1, optimise_fcn=fcn_list(2).name; end

% ========================== STRUCT = OPTIM.ALGO =========================
% optimisation banner
optim_report(optim_param.sys.output,' ');
optim_report(optim_param.sys.output,'===========================================');
optim_report(optim_param.sys.output,'=                                         =');
optim_report(optim_param.sys.output,'=               OPTIMISATION              =');
optim_report(optim_param.sys.output,'=                                         =');
optim_report(optim_param.sys.output,'===========================================');
optim_report(optim_param.sys.output,' ');

if isfield(optim,'extremum') && ismember(optim.extremum,{'maximum','minimum'}) && ~ismember(optim.method,{'grad_descent','grad_ascent'}) && ismember(optimise_fcn,{'fminnewton','fminkrotov'})
    optim_param.algo.extremum=optim.extremum;
    optim=rmfield(optim,'extremum'); % parsed -> remove
    optim_report(optim_param.sys.output,[pad('Find extremum',80) pad(optim_param.algo.extremum,20) ' (user-specified)']);
elseif ismember(optim.method,{'grad_descent'}) && ismember(optimise_fcn,{'fminnewton','fminkrotov'})
    optim_param.algo.extremum='minimum';
    optim_report(optim_param.sys.output,[pad('Find extremum',80) pad(optim_param.algo.extremum,20) ' (safe default)']);
elseif ismember(optim.method,{'grad_ascent'}) && ismember(optimise_fcn,{'fminnewton','fminkrotov'})
    optim_param.algo.extremum='maximum';
    optim_report(optim_param.sys.output,[pad('Find extremum',80) pad(optim_param.algo.extremum,20) ' (safe default)']);
elseif ismember(optimise_fcn,{'fminkrotov'})
    optim_param.algo.extremum='maximum';
    optim_report(optim_param.sys.output,[pad('Find extremum',80) pad(optim_param.algo.extremum,20) ' (safe default)']);
else
    optim_param.algo.extremum='minimum';
    optim_report(optim_param.sys.output,[pad('Find extremum',80) pad(optim_param.algo.extremum,20) ' (safe default)']);
end
switch optim_param.algo.extremum
    case 'maximum'
        optim_param.algo.max_min=+1;
    case 'minimum'
        optim_param.algo.max_min=-1;
end

if isfield(optim,'npen') && mod(optim.npen,1)==0
    optim_param.algo.pen=true;
    optim_param.algo.npen=optim.npen;
    optim=rmfield(optim,'npen'); % parsed -> remove    
    optim_report(optim_param.sys.output,[pad('Display penalty terms',80) pad('true',20) ' (user-specified)']);
    optim_report(optim_param.sys.output,[pad('Number of penalty terms',80) pad(int2str(optim_param.algo.npen),20) ' (user-specified)']);
else
    optim_param.algo.pen=false;
    optim_param.algo.npen=0;
end

% Optimisation method
switch optimise_fcn
    case {'fminnewton','fminkrotov'}
        if (~isfield(optim,'method') || ~ismember(optim.method,optim_algos)) && ~ismember(optimise_fcn,{'fminkrotov'})
            optim_param.algo.method='lbfgs';
            optim_report(optim_param.sys.output,[pad('Optimisation method',80) pad(optim_param.algo.method,20) ' (safe default)']);
        elseif (n_vars > lbfgs_switch) && ~ismember(optim.method,{'grad_descent','grad_ascent'}) && ~ismember(optimise_fcn,{'fminkrotov'})
            optim_param.algo.method='lbfgs';
            optim_report(optim_param.sys.output,[pad(['Optimisation method (variables > ' num2str(lbfgs_switch) ')'],80) pad(optim_param.algo.method,20) ' (safe default)']);
        elseif (~isfield(optim,'method') || ~ismember(optim.method,optim_algos))
            optim_param.algo.method='krotov';
            optim_report(optim_param.sys.output,[pad('Optimisation method',80) pad(optim_param.algo.method,20) ' (safe default)']);
        elseif ismember(optim.method,optim_algos) % must be stated explicitly for krotov
            optim_param.algo.method=optim.method;
            optim=rmfield(optim,'method'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Optimisation method',80) pad(optim_param.algo.method,20) ' (user-specified)']);
        end
        if ~ismember(optimise_fcn,{'fminkrotov'})
            % Store for the l-BFGS algorithm
            if ismember(optim_param.algo.method,{'lbfgs'})
                if isfield(optim,'n_store') && (mod(optim.n_store,1)==0) && optim.n_store>0 % if exists and is positive integer
                    optim_param.algo.n_store=optim.n_store; flag=0;
                    optim=rmfield(optim,'n_store'); % parsed -> remove
                else
                    optim_param.algo.n_store=20; flag=1;
                end
                % Trial if StoreN fits into memory
                succes=false;
                while(~succes)
                    try
                        optim_param.algo.deltaX=zeros(n_vars,optim_param.algo.n_store);
                        optim_param.algo.deltaG=zeros(n_vars,optim_param.algo.n_store);
                        succes=true;
                    catch ME
                        optim_report(optim_param.sys.output,[pad('WARNING: lbfgs store, out of memory...',80) pad(num2str(optim_param.algo.n_store),20) ' (decrease size)']);
                        succes=false;
                        optim_param.algo.deltaX=[]; optim_param.algo.deltaG=[];
                        optim_param.algo.n_store=optim_param.algo.n_store-1;
                        if(optim_param.algo.n_store<1)
                            rethrow(ME);
                        end
                    end
                end
                if flag, optim_report(optim_param.sys.output,[pad('lbfgs iteration store',80) pad(num2str(optim_param.algo.n_store),20) ' (safe default)']);
                else optim_report(optim_param.sys.output,[pad('lbfgs iteration store',80) pad(num2str(optim_param.algo.n_store),20) ' (user-specified)']);
                end
            end
        end
        
        if ismember(optim_param.algo.method,{'sr1','bfgs'}) && ~ismember(optimise_fcn,{'fminkrotov'})
            if isfield(optim,'inverse_method')
                optim_param.algo.inv_method=optim.inverse_method;
                optim=rmfield(optim,'inverse_method'); % parsed -> remove
                if optim_param.algo.inv_method, inverse_str='true'; else inverse_str='false'; end
                optim_report(optim_param.sys.output,[pad(['Inverse-' optim_param.algo.method ' Hessian update method'],80) pad(num2str(inverse_str),20) ' (user-specified)']);
            else
                optim_param.algo.inv_method=true;
                optim_report(optim_param.sys.output,[pad(['Inverse-' optim_param.algo.method ' Hessian update method'],80) pad('true',20) ' (safe default)']);
            end
        elseif ismember(optim_param.algo.method,{'sr1','bfgs','newton','krotov-bfgs','krotov-sr1'})
            optim_param.algo.inv_method=false;
        end
        
        % ========================== STRUCT = OPTIM.TOLS =========================
        % Maximum optimisation iterations
        if isfield(optim,'max_iterations') && (mod(optim.max_iterations,1)==0 || isinf(optim.max_iterations)) && optim.max_iterations>=0
            optim_param.tols.max_iter=optim.max_iterations;
            optim=rmfield(optim,'max_iterations'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Maximum optimisation iterations',80) pad(num2str(optim_param.tols.max_iter),20) ' (user-specified)']);
        else
            optim_param.tols.max_iter=100;
            optim_report(optim_param.sys.output,[pad('Maximum optimisation iterations',80) pad(num2str(optim_param.tols.max_iter),20) ' (safe default)']);
        end
        
        % termination tolerance on x
        if isfield(optim,'tol_x') && isnumeric(optim.tol_x)
            optim_param.tols.x=optim.tol_x;
            optim=rmfield(optim,'tol_x'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Termination tolerance on x',80) pad(num2str(optim_param.tols.x,'%0.8g'),20) ' (user-specified)']);
        else
            optim_param.tols.x=1e-6;
            optim_report(optim_param.sys.output,[pad('Termination tolerance on x',80) pad(num2str(optim_param.tols.x,'%0.8g'),20) ' (safe default)']);
        end
        
        if isfield(optim,'cost_display')
            optim_param.algo.f_disp=optim.cost_display;
            optim=rmfield(optim,'cost_display'); % parsed -> remove
        end
        
        % termination tolerance on f(x)
        if isfield(optim,'tol_fx') && isnumeric(optim.tol_fx)
            optim_param.tols.fx=optim.tol_fx;
            optim=rmfield(optim,'tol_fx'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Termination tolerance on f(x)',80) pad(num2str(optim_param.tols.fx,'%0.9g'),20) ' (user-specified)']);
        else
            optim_param.tols.fx=optim_param.algo.max_min*Inf;
        end
        
        % termination tolerance on the norm of gradient
        if isfield(optim,'tol_gfx') && isnumeric(optim.tol_gfx) && ~ismember(optimise_fcn,{'fminkrotov'})
            optim_param.tols.gfx=optim.tol_gfx;
            optim=rmfield(optim,'tol_gfx'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Termination tolerance on norm[grad(x)] (1st order optimility)',80) pad(num2str(optim_param.tols.gfx,'%0.8g'),20) ' (user-specified)']);
        elseif ~ismember(optimise_fcn,{'fminkrotov'})
            optim_param.tols.gfx=1e-6;
            optim_report(optim_param.sys.output,[pad('Termination tolerance on norm[grad(x)] (1st order optimility)',80) pad(num2str(optim_param.tols.gfx,'%0.8g'),20) ' (safe default)']);
        end
        
        % ======================= STRUCT = OPTIM.LINESEARCH ======================
        % linesearch does not apply to krotov methods
        if ~ismember(optimise_fcn,{'fminkrotov'})
            % add maximising/minimising strategy to linesreach structure
            switch optim_param.algo.extremum
                case 'maximum'
                    optim_param.linesearch.max_min=+1;
                case 'minimum'
                    optim_param.linesearch.max_min=-1;
            end
            
            % linesearch method
            if isfield(optim,'linesearch') && ismember(optim.linesearch,search_algos)
                optim_param.linesearch.method=optim.linesearch;
                optim=rmfield(optim,'linesearch'); % parsed -> remove
                optim_report(optim_param.sys.output,[pad('Line-search method',80) pad(optim_param.linesearch.method,20) ' (user-specified)']);
            else
                optim_param.linesearch.method='bracket-section';
                optim_report(optim_param.sys.output,[pad('Line-search method',80) pad(optim_param.linesearch.method,20) ' (safe default)']);
            end
            
            % no line-search is equivalent to 'newton-step'
            if ~ismember(optim_param.linesearch.method,{'newton-step'})
                % linesearch conditions
                if isfield(optim,'linesearch_rules') && ismember(optim.linesearch_rules,search_rules)
                    if ismember(optim.linesearch_rules,'Wolfe'), optim.linesearch_rules='Wolfe-strong'; end % Wolfe-strong is assumed for Wolfe
                    optim_param.linesearch.rules=optim.linesearch_rules;
                    optim=rmfield(optim,'linesearch_rules'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('Line-search conditions',80) pad(optim_param.linesearch.rules,20) ' (user-specified)']);
                elseif ismember(optim_param.linesearch.method,{'backtracking'})
                    optim_param.linesearch.rules='Armijo';
                    optim_report(optim_param.sys.output,[pad('Line-search conditions',80) pad(optim_param.linesearch.rules,20) ' (safe default)']);
                else
                    optim_param.linesearch.rules='Wolfe-strong';
                    optim_report(optim_param.sys.output,[pad('Line-search conditions',80) pad(optim_param.linesearch.rules,20) ' (safe default)']);
                end
                
                % Armijo condition: sufficient decrease
                if isfield(optim,'tol_linesearch_fx') && ismember(optim_param.linesearch.rules,{'Wolfe-strong','Wolfe-weak','Armijo'}) &&...
                        0<(optim.tol_linesearch_fx<1)
                    optim_param.linesearch.tols.c1=optim.tol_linesearch_fx;
                    optim=rmfield(optim,'tol_linesearch_fx'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('Armijo-Goldstein inequality: condition of sufficient decrease',80) pad(num2str(optim_param.linesearch.tols.c1,'%0.8g'),20) ' (user-specified)']);
                elseif isfield(optim,'tol_linesearch_fx') && ismember(optim_param.linesearch.rules,{'Goldstein'}) &&...
                        0<(optim.tol_linesearch_fx<0.5)
                    optim_param.linesearch.tols.c1=optim.tol_linesearch_fx;
                    optim=rmfield(optim,'tol_linesearch_fx'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('Armijo-Goldstein inequality: condition of sufficient decrease',80) pad(num2str(optim_param.linesearch.tols.c1,'%0.8g'),20) ' (user-specified)']);
                elseif ismember(optim_param.linesearch.rules,{'Wolfe-strong','Wolfe-weak','Armijo'}) && ismember(optim_param.algo.method,{'newton'})
                    optim_param.linesearch.tols.c1=1e-2;
                    optim_report(optim_param.sys.output,[pad('Armijo-Goldstein inequality: condition of sufficient decrease',80) pad(num2str(optim_param.linesearch.tols.c1,'%0.8g'),20) ' (safe default)']);
                elseif ismember(optim_param.linesearch.rules,{'Wolfe-strong','Wolfe-weak','Armijo'})
                    optim_param.linesearch.tols.c1=1e-2;
                    optim_report(optim_param.sys.output,[pad('Armijo-Goldstein inequality: condition of sufficient decrease',80) pad(num2str(optim_param.linesearch.tols.c1,'%0.8g'),20) ' (safe default)']);
                elseif ismember(optim_param.linesearch.rules,{'Goldstein'})
                    optim_param.linesearch.tols.c1=0.25;
                    optim_report(optim_param.sys.output,[pad('Armijo-Goldstein inequality: condition of sufficient decrease',80) pad(num2str(optim_param.linesearch.tols.c1,'%0.8g'),20) ' (safe default)']);
                end
                
                % Wolfe condition: curvature condition
                if isfield(optim,'tol_linesearch_gfx') && ismember(optim_param.linesearch.rules,{'Wolfe','Wolfe-strong'}) &&...
                        optim_param.linesearch.tols.c1<(optim.tol_linesearch_gfx<1)
                    optim_param.linesearch.tols.c2=optim.tol_linesearch_gfx;
                    optim=rmfield(optim,'tol_linesearch_gfx'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('Wolfe-Powel condition: curvature condition',80) pad(num2str(optim_param.linesearch.tols.c2,'%0.8g'),20) ' (user-specified)']);
                elseif ismember(optim_param.linesearch.rules,{'Wolfe-weak','Wolfe-strong'})
                    optim_param.linesearch.tols.c2=0.9;
                    optim_report(optim_param.sys.output,[pad('Wolfe-Powel condition: curvature condition',80) pad(num2str(optim_param.linesearch.tols.c2,'%0.8g'),20) ' (safe default)']);
                end
                
                % Backtracting linesearch contraction factor
                if isfield(optim,'tol_linesearch_reduction')&& ismember(optim_param.linesearch.method,{'backtracking'}) &&...
                        0<(optim.tol_linesearch_reduction<1)
                    optim_param.linesearch.tols.tau=optim.tol_linesearch_reduction;
                    optim=rmfield(optim,'tol_linesearch_reduction'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('Backtracting step reduction factor',80) pad(num2str(optim_param.linesearch.tols.tau,'%0.8g'),20) ' (user-specified)']);
                elseif ismember(optim_param.linesearch.method,{'backtracking'}) && ismember(optim_param.linesearch.rules,{'Goldstein'})
                    optim_param.linesearch.tols.tau=(2/(1+sqrt(5)));    %1/sqrt(2)
                    optim_report(optim_param.sys.output,[pad('Backtracting step reduction factor',80) pad(num2str(optim_param.linesearch.tols.tau,'%0.8g'),20) ' (safe default)']);
                elseif ismember(optim_param.linesearch.method,{'backtracking'})
                    optim_param.linesearch.tols.tau=(2/(1+sqrt(5)));    %0.5;%1/sqrt(2);
                    optim_report(optim_param.sys.output,[pad('Backtracting step reduction factor',80) pad(num2str(optim_param.linesearch.tols.tau,'%0.8g'),20) ' (safe default)']);
                end
                
                if ismember(optim_param.linesearch.method,{'bracket-section'})
                    % Bracket expansion if stepsize becomes larger
                    if isfield(optim,'tau1') && optim.tau1>1
                        optim_param.linesearch.tols.tau1=optim.tau1;
                        optim=rmfield(optim,'tau1'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad('Bracket expansion if stepsize becomes larger',80) pad(num2str(optim_param.linesearch.tols.tau1,'%0.8g'),20) ' (user-specified)']);
                    else
                        optim_param.linesearch.tols.tau1=3;
                        optim_report(optim_param.sys.output,[pad('Bracket expansion if stepsize becomes larger',80) pad(num2str(optim_param.linesearch.tols.tau1,'%0.8g'),20) ' (safe default)']);
                    end
                    
                    % Left bracket reduction used in section phase
                    if isfield(optim,'tau2') && 0<(optim.tau2<1)
                        optim_param.linesearch.tols.tau2=optim.tau2;
                        optim=rmfield(optim,'tau2'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad('Left bracket reduction used in section phase',80) pad(num2str(optim_param.linesearch.tols.tau2,'%0.8g'),20) ' (user-specified)']);
                    else
                        optim_param.linesearch.tols.tau2=0.1;
                        optim_report(optim_param.sys.output,[pad('Left bracket reduction used in section phase',80) pad(num2str(optim_param.linesearch.tols.tau2,'%0.8g'),20) ' (safe default)']);
                    end
                    
                    % Right bracket reduction used in section phase
                    if isfield(optim,'tau3') && 0<(optim.tau3<1)
                        optim_param.linesearch.tols.tau3=optim.tau3;
                        optim=rmfield(optim,'tau3'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad('Right bracket reduction used in section phase',80) pad(num2str(optim_param.linesearch.tols.tau3,'%0.8g'),20) ' (user-specified)']);
                    else
                        optim_param.linesearch.tols.tau3=0.5;
                        optim_report(optim_param.sys.output,[pad('Right bracket reduction used in section phase',80) pad(num2str(optim_param.linesearch.tols.tau3,'%0.8g'),20) ' (safe default)']);
                    end
                end
            end
        end
        
        % =========================== STRUCT = OPTIM.REG =========================
        % if newton, sr1, or bfgs method defined, set hessian tolerances
        if ismember(optim_param.algo.method,{'newton','sr1','bfgs','krotov-bfgs','krotov-sr1'})
            
            switch optim_param.algo.extremum
                case 'maximum'
                    optim_param.reg.def=-1;
                case 'minimum'
                    optim_param.reg.def=+1;
            end
            % Hessian regularisation method
            if isfield(optim,'regularisation') && ismember(optim.regularisation,reg_method)
                optim_param.reg.method=optim.regularisation;
                optim=rmfield(optim,'regularisation'); % parsed -> remove
                optim_report(optim_param.sys.output,[pad('Hessian/approximation regularisation method',80) pad(optim_param.reg.method,20) ' (user-specified)']);
            elseif ismember(optim_param.algo.method,{'newton'})
                optim_param.reg.method='RFO';
                optim_report(optim_param.sys.output,[pad('Hessian/approximation regularisation method',80) pad(optim_param.reg.method,20) ' (safe default)']);
            elseif ismember(optim_param.algo.method,{'sr1','bfgs'})
                optim_param.reg.method='';
                optim_report(optim_param.sys.output,[pad('Hessian/approximation regularisation method',80) pad('none',20) ' (safe default)']);
            elseif ismember(optim_param.algo.method,{'krotov-bfgs','krotov-sr1'})
                optim_param.reg.method='TRM';
                optim_report(optim_param.sys.output,[pad('Hessian/approximation regularisation method',80) pad(optim_param.reg.method,20) ' (safe default)']);
            else
                optim_param.reg.method='';
            end
            
            % Hessian conditioning method -
            if isfield(optim,'conditioning') && ismember(optim.conditioning,reg_cond)
                if ismember(optim_param.reg.method,{'RFO','TRM'})
                    if ismember(optim.conditioning,{'none',''})
                        optim_param.reg.cond_method='';
                        optim=rmfield(optim,'conditioning'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad('none',20) ' (user-specified)']);
                    else
                        optim_param.reg.cond_method=optim.conditioning;
                        optim=rmfield(optim,'conditioning'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad(optim_param.reg.cond_method,20) ' (user-specified)']);
                    end
                elseif ismember(optim_param.reg.method,{'CHOL'})
                    optim_param.reg.cond_method='';
                end
            elseif ismember(optim_param.reg.method,{'RFO'}) && any(ismember({'krotov'},optim_param.sys.costfun.fcn_calls))
                optim_param.reg.cond_method='scaled';
                optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad(optim_param.reg.cond_method,20) ' (safe default)']);
            elseif ismember(optim_param.reg.method,{'RFO'})
                optim_param.reg.cond_method='iterative';
                optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad(optim_param.reg.cond_method,20) ' (safe default)']);
            elseif ismember(optim_param.reg.method,{'TRM'})
                if ~isfield(optim,'max_cond'), optim_param.reg.cond_method='';
                else optim_param.reg.cond_method='iterative'; end
                optim_report(optim_param.sys.output,[pad('Hessian/approximation conditioning method',80) pad('none',20) ' (safe default)']);
            elseif ~ismember(optim_param.reg.method,{'CHOL'}) && ismember(optimise_fcn,{'fminkrotov'})
                optim_param.reg.method='RFO';
                optim_param.reg.cond_method='scaled';
                optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad(optim_param.reg.cond_method,20) ' (safe default)']);
            elseif ~ismember(optim_param.reg.method,{'CHOL'})
                optim_param.reg.method='RFO';
                optim_param.reg.cond_method='iterative';
                optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad(optim_param.reg.cond_method,20) ' (safe default)']);
            elseif isfield(optim,'max_cond')
                optim_param.reg.method='RFO';
                optim_param.reg.cond_method='iterative';
                optim_report(optim_param.sys.output,[pad(['Hessian/approximation conditioning method (' optim_param.reg.method ')'],80) pad(optim_param.reg.cond_method,20) ' (safe default)']);
                
            else
                optim_param.reg.cond_method='';
            end
            
            % Maximum Regularisation/conditioning iterates for RFO, TRM and CHOL
            if isfield(optim,'n_reg') && ismember(optim_param.reg.method,reg_method) && (mod(optim.n_reg,1)==0) && optim.n_reg>0
                optim_param.reg.n_reg=optim.n_reg;
                optim=rmfield(optim,'n_reg'); % parsed -> remove
                if ismember(optim_param.reg.method,{'CHOL'})
                    optim_report(optim_param.sys.output,[pad(['Maximum regularisation (' optim_param.reg.method ') iteratations'],80) pad(num2str(optim_param.reg.n_reg),20) ' (user-specified)']);
                elseif ismember(optim_param.reg.method,{'RFO','TRM'})
                    optim_report(optim_param.sys.output,[pad(['Maximum conditioning (' optim_param.reg.method ') iteratations'],80) pad(num2str(optim_param.reg.n_reg),20) ' (user-specified)']);
                end
            elseif ismember(optim_param.reg.method,reg_method)
                optim_param.reg.n_reg=2500;
                if ismember(optim_param.reg.method,{'CHOL'})
                    optim_report(optim_param.sys.output,[pad(['Maximum regularisation (' optim_param.reg.method ') iteratations'],80) pad(num2str(optim_param.reg.n_reg),20) ' (safe default)']);
                elseif ismember(optim_param.reg.method,{'RFO','TRM'})
                    optim_report(optim_param.sys.output,[pad(['Maximum conditioning (' optim_param.reg.method ') iteratations'],80) pad(num2str(optim_param.reg.n_reg),20) ' (safe default)']);
                end
            end
            
            if ~ismember(optim_param.reg.method,{'CHOL'})
                
                % Hessian conditioning method for RFO
                if isfield(optim,'alpha') && ismember(optim_param.reg.method ,{'RFO'}) && (optim.alpha>0)
                    optim_param.reg.alpha=optim.alpha;
                    optim=rmfield(optim,'alpha'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('RFO uniform scaling factor (alpha)',80) pad(num2str(optim_param.reg.alpha,'%0.8g'),20) ' (user-specified)']);
                elseif ismember(optim_param.reg.method,{'RFO'})
                    optim_param.reg.alpha=1;
                    optim_report(optim_param.sys.output,[pad('RFO uniform scaling factor (alpha)',80) pad(num2str(optim_param.reg.alpha,'%0.8g'),20) ' (safe default)']);
                end
                
                % Hessian delta value for TRM
                if isfield(optim,'delta') && ismember(optim_param.reg.method,{'TRM'}) && (optim.delta>0)
                    optim_param.reg.delta=optim.delta;
                    optim=rmfield(optim,'delta'); % parsed -> remove
                    optim_report(optim_param.sys.output,[pad('TRM uniform shifting factor (delta)',80) pad(num2str(optim_param.reg.delta,'%0.8g'),20) ' (user-specified)']);
                elseif ismember(optim_param.reg.method,{'TRM'})
                    optim_param.reg.delta=1;
                    optim_report(optim_param.sys.output,[pad('TRM uniform shifting factor (delta)',80) pad(num2str(optim_param.reg.delta,'%0.8g'),20) ' (safe default)']);
                end
                
                if ~isempty(optim_param.reg.cond_method)
                    % Hessian phi, alpha reduction factor for RFO
                    if isfield(optim,'phi') && ismember(optim_param.reg.method,{'RFO','TRM'}) && (optim.phi>0)
                        optim_param.reg.phi=optim.phi;
                        optim=rmfield(optim,'phi'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad([optim_param.reg.method ' conditioning factor (phi)'],80) pad(num2str(optim_param.reg.phi,'%0.8g'),20) ' (user-specified)']);
                    elseif ismember(optim_param.reg.method,{'RFO'})
                        optim_param.reg.phi=0.9;
                        optim_report(optim_param.sys.output,[pad([optim_param.reg.method ' conditioning factor (phi)'],80) pad(num2str(optim_param.reg.phi,'%0.8g'),20) ' (safe default)']);
                    elseif ismember(optim_param.reg.method,{'TRM'})
                        optim_param.reg.phi=1/(0.9)^2;
                        optim_report(optim_param.sys.output,[pad([optim_param.reg.method ' conditioning factor (phi)'],80) pad(num2str(optim_param.reg.phi,'%0.8g'),20) ' (safe default)']);
                    end
                    
                    % maximum condition number of Hessian for iterative or scaled conditioning
                    if isfield(optim,'max_cond') && ismember(optim_param.reg.method,{'RFO','TRM'}) && optim.max_cond>0
                        optim_param.reg.max_cond=optim.max_cond;
                        optim=rmfield(optim,'max_cond'); % parsed -> remove
                        optim_report(optim_param.sys.output,[pad('Maximum Hessian/approximation condition number',80) pad(num2str(optim_param.reg.max_cond,'%0.8e'),20) ' (user-specified)']);
                    elseif ismember(optim_param.reg.cond_method,{'iterative'})&&ismember(optim_param.reg.method,{'RFO','TRM'})
                        optim_param.reg.max_cond=1e4;
                        optim_report(optim_param.sys.output,[pad('Maximum Hessian/approximation condition number',80) pad(num2str(optim_param.reg.max_cond,'%0.8e'),20) ' (safe default)']);
                    elseif ismember(optim_param.reg.cond_method,{'scaled'})&&ismember(optim_param.reg.method,{'RFO','TRM'})
                        optim_param.reg.max_cond=1e14;
                        optim_report(optim_param.sys.output,[pad('Maximum Hessian/approximation condition number',80) pad(num2str(optim_param.reg.max_cond,'%0.8e'),20) ' (safe default)']);
                    end
                end
            end
            % optional eigenvalue tracking (should be used for investigating a problem or debugging)
            
        else
            optim_param.reg.method='';
            optim_param.reg.cond_method='';
        end
        
        if isfield(optim,'track_eig') && ismember(optim_param.algo.method,{'newton','sr1','bfgs','krotov-bfgs','krotov-sr1'})
            optim_param.reg.track_eig=optim.track_eig;
            optim=rmfield(optim,'track_eig'); % parsed -> remove
            if optim_param.reg.track_eig, track_str='true'; else track_str='false'; end
            optim_report(optim_param.sys.output,[pad('Explicitly track eigenvalues of Hessian/approximation',80) pad(num2str(track_str),20) ' (user-specified)']);
        else
            optim_param.reg.track_eig=false;
        end
        
    case {'fminsimplex'}
        
        if (~isfield(optim,'method') || ~ismember(optim.method,{'nm_simplex','md_simplex'}))
            optim_param.algo.method='nm_simplex';
            optim_report(optim_param.sys.output,[pad('Optimisation method',80) pad(optim_param.algo.method,20) ' (safe default)']);
        elseif ismember(optim.method,{'nm_simplex','md_simplex'}) % must be stated explicitly for krotov
            optim_param.algo.method=optim.method;
            optim=rmfield(optim,'method'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Optimisation method',80) pad(optim_param.algo.method,20) ' (user-specified)']);
        end
        
        if isfield(optim,'parallel_cost_fcn')
            optim_param.algo.parallel_cost_fcn=optim.parallel_cost_fcn;
            optim=rmfield(optim,'parallel_cost_fcn'); % parsed -> remove
            if optim_param.algo.parallel_cost_fcn, par_str='true'; else par_str='false'; end
            optim_report(optim_param.sys.output,[pad(['Parallel-' optim_param.algo.method ' cost function'],80) pad(num2str(par_str),20) ' (user-specified)']);
        else
            optim_param.algo.parallel_cost_fcn=false; % do not report to user - this is an advanced method
        end
        
        if (~isfield(optim,'simplex_min'))
            optim_param.tols.simplex_min=1e-3;
            optim_report(optim_param.sys.output,[pad('Tolerance for cgce test based on relative size of simplex',80) pad(num2str(optim_param.tols.simplex_min,'%0.8e'),20) ' (safe default)']);
        else
            optim_param.tols.simplex_min=optim.simplex_min;
            optim=rmfield(optim,'simplex_min'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Tolerance for cgce test based on relative size of simplex',80) pad(num2str(optim_param.tols.simplex_min,'%0.8e'),20) ' (user-specified)']);
        end
        % Maximum optimisation iterations
        if isfield(optim,'max_iterations') && (mod(optim.max_iterations,1)==0 || isinf(optim.max_iterations)) && optim.max_iterations>=0
            optim_param.tols.max_iter=optim.max_iterations;
            optim=rmfield(optim,'max_iterations'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Maximum optimisation iterations',80) pad(num2str(optim_param.tols.max_iter),20) ' (user-specified)']);
        else
            optim_param.tols.max_iter=100;
            optim_report(optim_param.sys.output,[pad('Maximum optimisation iterations',80) pad(num2str(optim_param.tols.max_iter),20) ' (safe default)']);
        end
        % maximum optimisation iterates
        if (~isfield(optim,'max_n_fx'))
            optim_param.tols.max_n_fx=Inf;
            optim_report(optim_param.sys.output,[pad('Max no. of f-evaluations',80) pad(num2str(optim_param.tols.max_n_fx),20) ' (safe default)']);
        else
            optim_param.tols.max_n_fx=optim.max_n_fx;
            optim=rmfield(optim,'max_n_fx'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Max no. of f-evaluations',80) pad(num2str(optim_param.tols.max_n_fx),20) ' (user-specified)']);
        end
        if (~isfield(optim,'termination'))
            optim_param.tols.termination=-Inf;
            optim_report(optim_param.sys.output,[pad('Target for f-values',80) pad(num2str(optim_param.tols.termination,'%0.8e'),20) ' (safe default)']);
        else
            optim_param.tols.termination=optim.termination;
            optim=rmfield(optim,'termination'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Target for f-values',80) pad(num2str(optim_param.tols.termination,'%0.8e'),20) ' (user-specified)']);
        end
        
        if (~isfield(optim,'init_simplex'))
            optim_param.algo.init_simplex='equilateral';
            optim_report(optim_param.sys.output,[pad('Initial simplex',80) pad(optim_param.algo.init_simplex,20) ' (safe default)']);
        else
            optim_param.algo.init_simplex=optim.init_simplex;
            optim=rmfield(optim,'init_simplex'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Initial simplex',80) pad(optim_param.algo.init_simplex,20) ' (user-specified)']);
        end
        
        if (~isfield(optim,'expansion'))
            optim_param.tols.expansion=2;
            optim_report(optim_param.sys.output,[pad('Expansion factor',80) pad(num2str(optim_param.tols.expansion,'%0.8e'),20) ' (safe default)']);
        else
            optim_param.tols.expansion=optim.expansion;
            optim=rmfield(optim,'expansion'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Expansion factor',80) pad(num2str(optim_param.tols.expansion,'%0.8e'),20) ' (user-specified)']);
        end
        if (~isfield(optim,'contraction'))
            optim_param.tols.contraction=0.5;
            optim_report(optim_param.sys.output,[pad('Contraction factor',80) pad(num2str(optim_param.tols.contraction,'%0.8e'),20) ' (safe default)']);
        else
            optim_param.tols.contraction=optim.contraction;
            optim=rmfield(optim,'contraction'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Contraction factor',80) pad(num2str(optim_param.tols.contraction,'%0.8e'),20) ' (user-specified)']);
        end
        
        if (~isfield(optim,'reflection')) && ismember(optim_param.algo.method,{'nm_simplex'})
            optim_param.tols.reflection=1.0;
            optim_report(optim_param.sys.output,[pad('Reflection factor',80) pad(num2str(optim_param.tols.reflection,'%0.8e'),20) ' (safe default)']);
        elseif ismember(optim_param.algo.method,{'nm_simplex'})
            optim_param.tols.reflection=optim.reflection;
            optim=rmfield(optim,'reflection'); % parsed -> remove
            optim_report(optim_param.sys.output,[pad('Reflection factor',80) pad(num2str(optim_param.tols.reflection,'%0.8e'),20) ' (user-specified)']);
        end
    otherwise
        error('unknown optimisation function')
end

% objective function display
if isfield(optim,'obj_fun_disp') && isa(optim.obj_fun_disp,'function_handle')
    obj_fun_disp_data=functions(optim.obj_fun_disp);
    if ismember(obj_fun_disp_data.type,{'anonymous'})
        optim_param.sys.obj_fun_disp=optim.obj_fun_disp;
        optim=rmfield(optim,'obj_fun_disp'); % parsed -> remove
        optim_report(optim_param.sys.output,[pad(['Objective function display = ', func2str(optim_param.sys.obj_fun_disp)],100) ' (user-specified)']);
        fx=sym('fx');
        obj_fun_str=char(optim_param.sys.obj_fun_disp(fx));
        obj_fun_str=strrep(obj_fun_str,'fx','f(x)');
        optim_param.sys.obj_fun_str=obj_fun_str;
    else
        error('objective function display should be an anonymous function')
    end
else
    % (default -> force)
    optim_param.sys.obj_fun_disp=@(f_x)(f_x);
    optim_param.sys.obj_fun_str='f(x)';
end
optim_param.sys.obj_fun_str=[optim_param.algo.extremum(1:3) '[' optim_param.sys.obj_fun_str ']'];

% Optimisation output function, evaluated at each iteration
if isa(optim_param.sys.outputfun,'function_handle')
    optim_report(optim_param.sys.output,[pad(['Optimisation output function = @', func2str(optim_param.sys.outputfun)],100) ' (user-specified)']);
end
unprocessed=fieldnames(optim);
if ~isempty(unprocessed)
    optim_report(optim_param.sys.output,'----------------------------------------------------------------------------------------------');
    for n=1:numel(unprocessed)
        optim_report(optim_param.sys.output,[pad('WARNING: unprocessed option',80) pad(unprocessed{n},20)]);
    end
    optim_report(optim_param.sys.output,'----------------------------------------------------------------------------------------------');
end



end

