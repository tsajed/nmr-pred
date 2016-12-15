% The waveforms on different channels are assumed to be stored in the 
% rows of the input array. The Hessian elements correspond to the ele-
% ments of the waveform array ordered as:
%
%                  [X1 Y1 Z1 X2 Y2 Z2 ... Xn Yn Zn]
%
% where X,Y,Z are different control channels and the index enumerates
% the time discretization points. Gradient dimensions and element or-
% der are the same as the input waveform dimensions and element order.
% Elements of the Hessian are reordered as to correspond to the wavef-
% orm array:
%                  [X1 X2 ... Xn Y1 Y2 ... Yn Z1 Z2 ... Zn] 
%
% interchanging the order from controls then time point to time point 
% then controls, or visa versa.
%
% Inputs:
%           hess    -   the old Hessian matrix to be reordered, curre-
%                       ntly ordered K first then N.
%
%           K       -   the first ordered variable of the old Hessian, 
%                       control channels in the example above.
%
%           N       -   the second ordered variable of the old Hessian,
%                       time points in the example above.
%
% Output:
%           hess_new -  the newly order Hessian with N first then K.
% 
% d.goodwin@soton.ac.uk

function hess_new=hess_reorder(hess,K,N)

% Check consistency
grumble(hess,K,N)

% Pre-allocate space for new Hessian
hess_new=zeros(size(hess));

% Loop over the elements
parfor j=1:(N*K)^2 
    
    % Dummy variable for parfor to work
    hess_temp=zeros(size(hess));
    
    % Row and column of new Hessian
    row=1+mod(j-1,N*K);
    col=1+(j-1-mod(j-1,N*K))/(N*K);
    
    % Variables corresponding to elements
    k1=1+mod(row-1,N); k2=1+mod(col-1,N);
    n1=1+(row-1-mod(row-1,N))/(N);
    n2=1+(col-1-mod(col-1,N))/(N); 
    
    % Row and column of the old Hessian
    row_old=n1+K*(k1-1);
    col_old=n2+K*(k2-1);
    
    % Allocate the old Hessian to its new place
    hess_temp(row,col)=hess(row_old,col_old);
    hess_new=hess_new+hess_temp;
    
end

end

% Consistency enforcement
function grumble(hess,dim1,dim2)
if numel(hess)~=(dim1*dim2)^2
    error('Hessian size should be (K*N)^2')
end
end

% If it's true that our species is alone in the universe, 
% then I'd have to say the universe aimed rather low and
% settled for very little.
%
% George Carlin

