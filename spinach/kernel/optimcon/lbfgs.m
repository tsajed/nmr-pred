% L-BFGS optimization step update.
%
% <http://spindynamics.org/wiki/index.php?title=Lbfgs.m>

function direction=lbfgs(x_hist,df_hist,grad,N)

% Initialize variables
a=zeros(1,N);
p=zeros(1,N);

for i=1:N
    p(i)= 1 / (df_hist(:,i)'*x_hist(:,i));
    a(i) = p(i)* x_hist(:,i)' * grad;
    grad = grad - a(i) * df_hist(:,i);
end
% Scaling of initial Hessian (identity matrix)
p_k = df_hist(:,1)'*x_hist(:,1) / sum(df_hist(:,1).^2);

% Make r = - Hessian * gradient
direction = p_k * grad;
for i=N:-1:1,
    b = p(i) * df_hist(:,i)' * direction;
    direction = direction + x_hist(:,i)*(a(i)-b);
end

direction=-direction;

end

% If someone steals your password, you can change it. But if someone
% steals your thumbprint, you can't get a new thumb. The failure mod-
% es are very different.
%
% Bruce Schneier, on biometric security

