% function to call and collect the correct amount of outputs from an 
% objective function
%
% <http://spindynamics.org/wiki/index.php?title=Objeval.m>

function [data,fx,grad,hess]=objeval(x,obj_fun_handle, data)
% Call the cost function with gradient and Hessian if nessesary
fcn_list=dbstack;
if numel(fcn_list)==1, optimise_fcn={}; end
if numel(fcn_list)>1, optimise_fcn=fcn_list(2).name; end

if ismember(optimise_fcn,{'linesearch','fminnewton','sectioning','bracketing'})
    timem=tic;
    if ( nargout ==2 )
        if ~isempty(data.cost_fun_vars)
            [fx,data.cost_fun_vars]=...
                feval(obj_fun_handle,reshape(x,data.x_shape),data.cost_fun_vars);
        else
            [fx]=feval(obj_fun_handle,reshape(x,data.x_shape));
        end
        data.time_fx=data.time_fx+toc(timem);
        data.n_fx=data.n_fx+1;
    elseif ( nargout ==3 )
        if ~isempty(data.cost_fun_vars)
            [fx,grad,data.cost_fun_vars]=...
                feval(obj_fun_handle,reshape(x,data.x_shape),data.cost_fun_vars);
        else
            [fx,grad]=feval(obj_fun_handle,reshape(x,data.x_shape));
        end
        data.time_gfx=data.time_gfx+toc(timem);
        data.n_fx=data.n_fx+1; data.n_grad=data.n_grad+1;
        grad=grad(:);
    elseif ( nargout ==4 )
        if ~isempty(data.cost_fun_vars)
            [fx,grad,hess,data.cost_fun_vars]=...
                feval(obj_fun_handle,reshape(x,data.x_shape),data.cost_fun_vars);
        else
            [fx,grad,hess]=feval(obj_fun_handle,reshape(x,data.x_shape));
        end
        data.time_hfx=data.time_hfx+toc(timem);
        data.n_fx=data.n_fx+1; data.n_grad=data.n_grad+1; data.n_hess=data.n_hess+1;
        grad=grad(:);
    end
    
    % store then sum fx if penalties are included
    if data.npen>0, data.fx_sep_pen=fx; fx=sum(fx); end
    
elseif ismember(optimise_fcn,{'fminsimplex'})
    
    fx=zeros(size(x,2),data.npen+1);
    if data.parallel_cost_fcn
        X=cell(size(x,2),1);
        for n=1:numel(X), X{n}=reshape(x(:,n),data.x_shape); end
        
        % Call the cost function with gradient and Hessian if nessesary
        timem=tic;
        if isempty(data.cost_fun_vars)
            [Fx]=feval(obj_fun_handle,X);
        else
            [Fx,data.cost_fun_vars]=feval(obj_fun_handle,X,data.cost_fun_vars);
        end
        for n=1:numel(X), fx(n,:)=Fx{n}; end
        data.time_fx=data.time_fx+toc(timem);
        data.n_fx=data.n_fx+numel(X);
    else
        for n=1:size(x,2)
            % Call the cost function with gradient and Hessian if nessesary
            timem=tic;
            if isempty(data.cost_fun_vars)
                [fx(n,:)]=feval(obj_fun_handle,reshape(x(:,n),data.x_shape));
            else
                [fx(n,:),data.cost_fun_vars]=feval(obj_fun_handle,reshape(x(:,n),data.x_shape),data.cost_fun_vars);
            end
            data.time_fx=data.time_fx+toc(timem);
            data.n_fx=data.n_fx+1;
        end
    end
    % store then sum fx if penalties are included
    if data.npen>0, data.fx_sep_pen=fx; fx=sum(fx,2); end
end
end