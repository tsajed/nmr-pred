% Performs a line search to find an appropriate step length called from 
% within a given optimisation algorithm, e.g. using fminnewton.m. This
% function generally assumes cheap gradients.
%
% bracket-section line search based on fminlbfgs.m code from D. Kroon,
% University of Twente (Nov 2010).
%
% <http://spindynamics.org/wiki/index.php?title=Linesearch.m>

function [alpha,fx_1,gfx_1,exitflag,data]=linesearch(cost_function,...
    d_0, x_0, fx_0, gfx_0, data, optim)

% set alpha to alpha_0
if (ismember(optim.algo.method,{'grad_ascent','grad_descent'})) ||...
        (ismember(optim.algo.method,{'lbfgs','bfgs','sr1'}) && data.iter==0)
    alpha=min(1/norm(gfx_0,Inf),5);
else
    alpha=1;
end

% set new x to that with initial alpha
x_1=x_0+alpha*d_0;

% preallocate initial functional and gradient
fx_1=[]; gfx_1=[]; exitflag=[];

switch optim.linesearch.method
    
    case {'bracket-section'}
        
        % initial evaluation of line search
        switch optim.linesearch.rules
            case {'Armijo','Goldstein'}
                [data,fx_1]=objeval(x_1,cost_function, data);
            case {'Wolfe','Wolfe-weak','Wolfe-strong'}
                [data,fx_1,gfx_1]=objeval(x_1,cost_function, data);
        end
        
        % Find a bracket [A B] of acceptable points
        [A, B, alpha, fx_1, gfx_1, exitflag, data] = bracketing(cost_function,...
            alpha, d_0, x_0, fx_0, gfx_0, fx_1, gfx_1, data, optim.linesearch);
        
        if (exitflag  == 2)
            
            % Find acceptable point within bracket
            [alpha, fx_1, gfx_1, exitflag, data] = sectioning(cost_function,...
                A, B, x_0, fx_0, gfx_0, d_0, data, optim.linesearch);
        end
        
        
    case {'backtracking'}
        
        fminimum = fx_0 - 1e16*(1+abs(fx_0));
        
        % find acceptable search length
        while (true)
            
            % initial evaluation of line search
            switch optim.linesearch.rules
                case {'Armijo','Goldstein'}
                    [data,fx_1]=objeval(x_1,cost_function, data);
                case {'Wolfe','Wolfe-weak','Wolfe-strong'}
                    [data,fx_1,gfx_1]=objeval(x_1,cost_function, data);
            end
            
            % Terminate if f < fminimum
            if (fx_1 <= fminimum), exitflag = 3; return; end
            
            % Goldstein conditions of sufficient decrease satisfied
            if alpha_conditions(1,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,optim.linesearch) &&...
                    alpha_conditions(2,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,optim.linesearch)
                break,
            end
            
            % set new alpha with contraction factor (rho)
            alpha = alpha*optim.linesearch.tols.tau;
            
            % check step length is not too small
            if alpha < eps, exitflag = -2; break, end
            
            % update x with new alpha
            x_1=x_0+alpha*d_0;
            
        end
        
    case {'newton-step'}
        
        %already done
        alpha=1; fx_1=[]; gfx_1=[];
        
    otherwise
        
        error('unknown line search method')
        
end

end

function [A, B, alpha, fx, gfx, exitflag, data] = bracketing(cost_function,...
    alpha, d_0, x_0, fx_0, gfx_0, fx_2, gfx_2, data, ls_param)

% Initialise bracket A and bracket B
A.alpha = []; A.fx = []; A.gfx = [];
B.alpha = []; B.fx = []; B.gfx = [];

% Set maximum value of alpha (determined by fminimum)
fminimum = (-1e16*(1+abs(fx_0)) - ls_param.max_min*fx_0);
alpha_max = ls_param.max_min*(fx_0 - fminimum)/(ls_param.tols.c1*gfx_0'*d_0);

% Evaluate f(alpha) and f'(alpha)
fx=fx_0;        fx_1 = fx_0;
gfx = gfx_0;	gfx_1 = gfx_0;
alpha_1=0;      alpha_2=alpha;

while(true)
    
    % Terminate if f < fminimum
    if (fx_2 <= fminimum), 
        exitflag = 3; return; 
    end
    
    % Bracket located - case 1 (Wolfe conditions)
    if ~alpha_conditions(1,alpha,fx_0,fx_2,gfx_0,[],d_0,ls_param) ||...
            ~alpha_conditions(0,[],fx_0,fx_2,[],[],[],ls_param)
        % Set the bracket values
        A.alpha = alpha_1; A.fx = fx_1;  A.gfx = gfx_1;
        B.alpha = alpha_2; B.fx = fx_2;  B.gfx = gfx_2;
        % Finished bracketing phase
        exitflag  = 2; return
    end
    
    % Acceptable steplength found
    if alpha_conditions(2,alpha,fx_0,fx_2,gfx_0,gfx_2,d_0,ls_param)
        
        % Store the found alpha values
        alpha=alpha_2; fx=fx_2; gfx=gfx_2;
        
        % Finished bracketing phase, and no need to call sectioning phase
        exitflag = [];  return
    end
    
    % Bracket located - case 2
    if ~alpha_conditions(3,[],[],[],[],gfx_2,d_0,ls_param)
        % Set the bracket values
        A.alpha = alpha_2; A.fx = fx_2;  A.gfx = gfx_2;
        B.alpha = alpha_1; B.fx = fx_1;  B.gfx = gfx_1;
        % Finished bracketing phase
        exitflag  = 2; return
    end
    
    % Update alpha (2*alpha_2 - alpha_1 > alpha_max )
    if (2*alpha_2 - alpha_1 < alpha_max )
        brcktEndpntA = 2*alpha_2-alpha_1;
        brcktEndpntB = min(alpha_max,alpha_2+ls_param.tols.tau1*(alpha_2-alpha_1));
        % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree polynomial
        % that interpolates f() and f'() at alphaPrev and at alpha
        alpha_new = cubic_interpolation(brcktEndpntA,brcktEndpntB,...
            alpha_1, alpha_2, fx_1, gfx_1'*d_0, fx_2, gfx_2'*d_0,ls_param.max_min);
        alpha_1 = alpha_2;
        alpha_2=alpha_new;
    else
        alpha_2 = alpha_max;
    end
    % Evaluate f(alpha) and f'(alpha)
    fx_1 = fx_2; gfx_1 = gfx_2;
    x_1=x_0+alpha_2*d_0;
    
    % Calculate value and gradient of current alpha
    [data,fx_2, gfx_2]=objeval(x_1,cost_function, data);
end

end

function [alpha, fx_1, gfx_1, exitflag, data] = sectioning(cost_function, A, B, x_0, fx_0, gfx_0, d_0, data, ls_param)

alpha=[]; fx_1=[]; gfx_1=[];
while(true)
    
    % Pick alpha in reduced bracket
    End_A = A.alpha + min(ls_param.tols.tau2,ls_param.tols.c2)*(B.alpha - A.alpha);
    End_B = B.alpha - ls_param.tols.tau3*(B.alpha - A.alpha);
    
    % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree
    % polynomial that interpolates f() and f'() at "a" and at "b".
    alpha = cubic_interpolation(End_A,End_B,...
        A.alpha, B.alpha, A.fx, A.gfx'*d_0, B.fx, B.gfx'*d_0,ls_param.max_min);
    
    % No acceptable point could be found
    if (abs( (alpha - A.alpha)*(A.gfx'*d_0) ) <= eps(max(1,abs(fx_0))))
        exitflag = -2; return;
    end
    
    % Calculate value and gradient of current alpha
    [data,fx_1, gfx_1]=objeval(x_0+alpha*d_0, cost_function, data);
    
    % Store current bracket position of A
    Tmp=A;
    
    % Update the current brackets
    if ~alpha_conditions(1,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,ls_param) ||...
            ~alpha_conditions(0,alpha,A.fx,fx_1,A.gfx,gfx_1,d_0,ls_param)
        
        % Update bracket B to current alpha
        B.alpha = alpha; B.fx = fx_1; B.gfx = gfx_1;
    else
        % Wolfe conditions, if true then acceptable point found
        if alpha_conditions(2,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,ls_param)
            
            exitflag = []; return,
        end
        
        % Update bracket A
        A.alpha = alpha; A.fx = fx_1;  A.gfx = gfx_1;
        
        % B becomes old bracket A;
        if ls_param.max_min*(A.alpha - B.alpha)*(gfx_1'*d_0) >= 0, B=Tmp; end
        
    end
    
    % No acceptable point could be found
    if (abs(B.alpha-A.alpha) < eps), exitflag = -2; return, end
    
end

end

function [alpha,fx]= cubic_interpolation(End_A,End_B,alpha_A,alpha_B,f_A,dir_deriv_A,f_B,dir_deriv_B,max_min)

% determines the coefficients of the cubic polynomial
c1=-2*(f_B-f_A) + (dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c2= 3*(f_B-f_A) - (2*dir_deriv_A+dir_deriv_B)*(alpha_B-alpha_A);
c3=(alpha_B-alpha_A)*dir_deriv_A;
c4=f_A;

% Convert bounds to the z-space
bounds = ([End_A End_B ]-alpha_A)./(alpha_B - alpha_A);

% Find minima and maxima from the roots of the derivative of the polynomial.
sPoints = roots([3*c1 2*c2 1*c3]);

% Remove imaginary and points outside range and make vector with solutions
sPoints(imag(sPoints)~=0)=[];
sPoints(sPoints<min(bounds))=[];
sPoints(sPoints>max(bounds))=[];
sPoints=[min(bounds) sPoints(:)' max(bounds)];

% Select the global minimum point
[fx,k]=max(max_min*polyval([c1 c2 c3 c4],sPoints)); fx=max_min*fx;

% Add the offset and scale back from [0..1] to the alpha domain
alpha = alpha_A + sPoints(k)*(alpha_B - alpha_A);

end

function test=alpha_conditions(test_type,alpha,fx_0,fx_1,gfx_0,gfx_1,d_0,ls_def)
% note: ls_def.max_min switches maximisation (+1) and minimisation (-1)
fx_0=ls_def.max_min.*fx_0;
fx_1=ls_def.max_min.*fx_1;

switch ls_def.rules
    case {'Armijo'}
        if test_type==0
            
            % increase (max) or decrease (min) of functional value -
            % gradient free test
            test=(fx_1 > fx_0);
            
        elseif test_type==1
            
            % sufficient increase( max) or decrease (min) of function -
            % also called Armijo condition
            test=(fx_1 >= fx_0 + ls_def.max_min.*ls_def.tols.c1*alpha*(gfx_0'*d_0));
            
        elseif test_type==2
            test=1;
        elseif test_type==3
            test=1;
        end
    case {'Goldstein'}
        if test_type==0
            
            % increase (max) or decrease (min) of functional value -
            % gradient free test
            test=(fx_1 > fx_0);
            
        elseif test_type==1
            
            % sufficient increase( max) or decrease (min) of function -
            % also called Armijo condition
            test=(fx_1 >= fx_0 + ls_def.max_min.*ls_def.tols.c1*alpha*(gfx_0'*d_0));
            
        elseif test_type==2
            test=(fx_1 >= fx_0 + ls_def.max_min.*(1-ls_def.tols.c1)*alpha*(gfx_0'*d_0));
        elseif test_type==3
            test=(gfx_1'*d_0 < 0);
        end
    case {'Wolfe-weak','Wolfe'}
        if test_type==0
            
            % increase (max) or decrease (min) of functional value -
            % gradient free test
            test=(fx_1 > fx_0);
            
        elseif test_type==1
            
            % sufficient increase( max) or decrease (min) of function -
            % also called Armijo condition
            test=(fx_1 >= fx_0 + ls_def.max_min.*ls_def.tols.c1*alpha*(gfx_0'*d_0));
            
        elseif test_type==2
            test=(gfx_1'*d_0 >= ls_def.tols.c2*(gfx_0'*d_0));
        elseif test_type==3
            test=(ls_def.max_min.*gfx_1'*d_0 > 0);
        end
    case {'Wolfe-strong'}
        if test_type==0
            
            % increase (max) or decrease (min) of functional value -
            % gradient free test
            test=(fx_1 > fx_0);
            
        elseif test_type==1
            
            % sufficient increase( max) or decrease (min) of function -
            % also called Armijo condition
            test=(fx_1 >= fx_0 + ls_def.max_min.*ls_def.tols.c1*alpha*(gfx_0'*d_0));
            
        elseif test_type==2
            test=(abs(gfx_1'*d_0) <= ls_def.tols.c2*abs(gfx_0'*d_0));
        elseif test_type==3
            test=(ls_def.max_min.*gfx_1'*d_0 > 0);
        end
    otherwise
        error('unknown linesearch conditions')
end

end

