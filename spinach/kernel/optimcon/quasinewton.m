% BFGS and SR1 Hessian update.
%
% <http://spindynamics.org/wiki/index.php?title=Quasinewton.m>

function B_new=quasinewton(B_old,x_new,x_old,df_new,df_old,inv_method,algo)

% Check consistency
grumble(x_new,x_old,df_new,df_old);

% if no method is supplied, use the bfgs formula
if ~exist('algo','var'), algo='bfgs'; end

% default inv_method is true if not given or not logical
if ~exist('inv_method','var'), inv_method=true; end
if ~islogical(inv_method), inv_method=true; end

% Eliminate rounding errors
B_old=real((B_old+B_old')/2);

% Compute shorthands
dx=x_new-x_old; ddf=df_new-df_old;

switch algo
    
    case {'bfgs'}
        
        if inv_method
            
            % Update the inverse Hessian
            B_new=B_old+(dx'*ddf+ddf'*B_old*ddf)*(dx*dx')/((dx'*ddf)^2)-...
                ((B_old*ddf)*dx'+dx*(ddf'*B_old))/(dx'*ddf);
            
        elseif ~inv_method
            
            % Update the Hessian
            B_new=B_old+(ddf*ddf')/(ddf'*dx)-((B_old*dx)*(dx'*B_old))/(dx'*B_old*dx);
            
        end
        
    case {'sr1'}
        
        if inv_method
            
            % Update the inverse Hessian
            B_new=B_old+((dx-B_old*ddf)*(dx-B_old*ddf)')/((dx-B_old*ddf)'*ddf);
            
        elseif ~inv_method
            
            % Update the Hessian
            B_new=B_old+((ddf-B_old*dx)*(ddf-B_old*dx)')/((ddf-B_old*dx)'*dx);
            
        end
        
    otherwise
        error('unsupported quasi-newton update method')
end

% Eliminate rounding errors
B_new=real((B_new+B_new')/2);

end

% Consistency enforcement
function grumble(x_new,x_old,df_new,df_old)
if (size(x_new,2)~=1)||(size(x_old,2)~=1)||(size(df_old,2)~=1)||(size(df_new,2)~=1)
    error('all vector inputs must be column vectors.');
end
if (~isreal(x_new))||(~isreal(x_old))||(~isreal(df_new))||(~isreal(df_old))
    error('all vector inputs must be real.');
end
end

% Of all tyrannies, a tyranny exercised for the good of its victims may be
% the most oppressive. It may be better to live under robber barons than
% under omnipotent moral busybodies. The robber baron's cruelty may some-
% times sleep, his cupidity may at some point be satiated; but those who
% torment us for our own good will torment us without end, for they do so
% with the approval of their consciences.
%
% C.S. Lewis

