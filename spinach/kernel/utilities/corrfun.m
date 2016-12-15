% Rotational correlation function, normalized to be the correlation func-
% tion between second rank Wigner functions. Syntax:
%
%              [weights,rates]=corrfun(spin_system,k,m,p,q)
%
% where the indices k,m,p,q correspond to the four indices found in the
% ensemble-averaged Wigner function product:
%
%                       G(k,m,p,q)=<D2(k,m)*D2(p,q)'>
% 
% the function requires spin_system.rlx.tau_c to be either a single cor-
% relation time (in which case the isotropic rotational diffusion model
% is used) or vector with two correlation times (in which case they are
% assumed to be correlation times for rotation around and perpendicular-
% ly to the main axis respectively) or vector with three correlation ti-
% mes (in which case those are assumed to be the correlation times for
% the rotation around the XX, YY and ZZ direction respectively of the
% rotational diffusion tensor).
%
% The function returns the weights and decay rates of the individual ex-
% ponentials in the correlation function. Note that the decay rates re-
% turned are negative numbers.
%
% Note: Wigner function indices are sorted in descending order, that is,
%       [1 2 3 4 5] in the input maps onto [2 1 0 -1 -2].
%
% i.kuprov@soton.ac.uk

function [weights,rates]=corrfun(spin_system,k,m,p,q)

% Check consistency
grumble(k,m,p,q);

% Select rotational diffusion model
switch numel(spin_system.rlx.tau_c)
    
    case 1
        
        % Use isotropic rotational diffusion model
        D=1/(6*spin_system.rlx.tau_c); rates=-6*D;
        weights=(1/5)*krondelta(k,p)*krondelta(m,q);
         
    case 2
        
        % Use axial rotational diffusion model
        D_ax=1/(6*spin_system.rlx.tau_c(1));
        D_eq=1/(6*spin_system.rlx.tau_c(2));
        weights=(1/5)*krondelta(k,p)*krondelta(m,q);
        rates=-(6*D_eq+((m-3)^2)*(D_ax-D_eq));
       
    case 3
        
        % Use anisotropic rotational diffusion model
        Dxx=1/(6*spin_system.rlx.tau_c(1));
        Dyy=1/(6*spin_system.rlx.tau_c(2));
        Dzz=1/(6*spin_system.rlx.tau_c(3));
        
        % Refuse to process degenerate cases
        if (abs(Dxx-Dyy)<1e-6*mean([Dxx Dyy Dzz]))||...
           (abs(Dyy-Dzz)<1e-6*mean([Dxx Dyy Dzz]))||...
           (abs(Dzz-Dxx)<1e-6*mean([Dxx Dyy Dzz]))
            error('the three rotational correlation times must be different.');
        end
        
        % Compute decay rates
        delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz);
        rates(1)=-(4*Dxx+Dyy+Dzz);
        rates(2)=-(Dxx+4*Dyy+Dzz);
        rates(3)=-(Dxx+Dyy+4*Dzz);
        rates(4)=-(2*Dxx+2*Dyy+2*Dzz-2*delta);
        rates(5)=-(2*Dxx+2*Dyy+2*Dzz+2*delta);
        
        % Compute coefficients
        lambda_p=sqrt(2/3)*(Dxx+Dyy-2*Dzz+2*delta)/(Dxx-Dyy);
        lambda_m=sqrt(2/3)*(Dxx+Dyy-2*Dzz-2*delta)/(Dxx-Dyy);
        h(1,2)=1/sqrt(2); h(1,4)=1/sqrt(2); h(2,2)=-1/sqrt(2);
        h(2,4)=1/sqrt(2); h(3,5)=1/sqrt(2); h(3,1)=-1/sqrt(2); 
        h(4,1)=1/sqrt(2+lambda_m^2); h(4,3)=lambda_m/sqrt(2+lambda_m^2); h(4,5)=1/sqrt(2+lambda_m^2);
        h(5,1)=1/sqrt(2+lambda_p^2); h(5,3)=lambda_p/sqrt(2+lambda_p^2); h(5,5)=1/sqrt(2+lambda_p^2);
        
        % Compute weights 
        for j=1:5, weights(j)=(1/5)*krondelta(k,p)*h(j,m)*h(j,q); end %#ok<AGROW>
        
end

end

% Consistency enforcement
function grumble(k,m,p,q)
if (~isnumeric(k))||(~isnumeric(m))||...
   (~isnumeric(p))||(~isnumeric(q))
    error('all indices must be numeric.');
end
if ~isreal([k m p q])
    error('all indices must be real.');
end
if any(mod([k m p q],1)~=0)
    error('all indices must be integers.');
end
if any([k m p q]>5)||any([k m p q]<1)
    error('k,m,p,q must be from [1,5] interval.');
end
end

% "But they are useless. They can only give you answers."
%
% Pablo Picasso, about computers.

