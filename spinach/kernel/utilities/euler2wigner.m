% Second-rank Wigner rotation matrix as a function of Euler angles. Two
% possible input styles are:
%
%                    W=euler2wigner(alpha,beta,gamma)
%                    W=euler2wigner([alpha beta gamma])
%
% where alpha, beta and gamma are Euler angles in radians. Rows and columns
% of the resulting Wigner matrix are sorted by descending ranks, that is:
%    
%                         [D( 2,2)  ...  D( 2,-2)
%                            ...    ...    ...  
%                          D(-2,2)  ...  D(-2,-2)]
%
% The resulting Wigner matrix is to be used as v=W*v, where v is a column
% vector of irreducible spherical tensor coefficients, listed vertically in
% the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).
% 
% i.kuprov@soton.ac.uk
% olga.efimova@chem.ox.ac.uk

function W=euler2wigner(arg1,arg2,arg3)

% Adapt to the input style
if nargin==1
    % Assume that a single input is a 3-vector
    alp=arg1(1); bet=arg1(2); gam=arg1(3);
elseif nargin==3
    % Assume that three inputs are Euler angles
    alp=arg1; bet=arg2; gam=arg3;
else
    % Bomb out in all other cases
    error('incorrect number of input arguments.');
end

% Check consistency
grumble(alp,bet,gam);

% Do the math, as per Brink and Satchler, Eq 2.17
W=[ exp(-1i*(+2)*alp)*  (+cos(bet/2)^4)*                   exp(-1i*(+2)*gam), ...
    exp(-1i*(+2)*alp)*  (-sin(bet)*(cos(bet)+1)/2)*        exp(-1i*(+1)*gam), ...
    exp(-1i*(+2)*alp)*  (+sqrt(3/8)*sin(bet)^2)*           exp(-1i*( 0)*gam), ...
    exp(-1i*(+2)*alp)*  (+sin(bet)*(cos(bet)-1)/2)*        exp(-1i*(-1)*gam), ...
    exp(-1i*(+2)*alp)*  (+sin(bet/2)^4)*                   exp(-1i*(-2)*gam); ...

    exp(-1i*(+1)*alp)*  (+sin(bet)*(cos(bet)+1)/2)*        exp(-1i*(+2)*gam), ...
    exp(-1i*(+1)*alp)*  (+(2*cos(bet)-1)*(cos(bet)+1)/2)*  exp(-1i*(+1)*gam), ...
    exp(-1i*(+1)*alp)*  (-sqrt(3/2)*sin(bet)*cos(bet))*    exp(-1i*( 0)*gam), ...
    exp(-1i*(+1)*alp)*  (-(2*cos(bet)+1)*(cos(bet)-1)/2)*  exp(-1i*(-1)*gam), ...
    exp(-1i*(+1)*alp)*  (+sin(bet)*(cos(bet)-1)/2)*        exp(-1i*(-2)*gam); ...

    exp(-1i*( 0)*alp)*  (+sqrt(3/8)*sin(bet)^2)*           exp(-1i*(+2)*gam), ...
    exp(-1i*( 0)*alp)*  (+sqrt(3/2)*sin(bet)*cos(bet))*    exp(-1i*(+1)*gam), ...
    exp(-1i*( 0)*alp)*  (3*cos(bet)^2-1)/2*                exp(-1i*( 0)*gam), ...
    exp(-1i*( 0)*alp)*  (-sqrt(3/2)*sin(bet)*cos(bet))*    exp(-1i*(-1)*gam), ...
    exp(-1i*( 0)*alp)*  (+sqrt(3/8)*sin(bet)^2)*           exp(-1i*(-2)*gam); ...

    exp(-1i*(-1)*alp)*  (-sin(bet)*(cos(bet)-1)/2)*        exp(-1i*(+2)*gam), ...
    exp(-1i*(-1)*alp)*  (-(2*cos(bet)+1)*(cos(bet)-1)/2)*  exp(-1i*(+1)*gam), ...
    exp(-1i*(-1)*alp)*  (+sqrt(3/2)*sin(bet)*cos(bet))*    exp(-1i*( 0)*gam), ...
    exp(-1i*(-1)*alp)*  (+(2*cos(bet)-1)*(cos(bet)+1)/2)*  exp(-1i*(-1)*gam), ...
    exp(-1i*(-1)*alp)*  (-sin(bet)*(cos(bet)+1)/2)*        exp(-1i*(-2)*gam); ...

    exp(-1i*(-2)*alp)*  (+sin(bet/2)^4)*                   exp(-1i*(+2)*gam), ...
    exp(-1i*(-2)*alp)*  (-sin(bet)*(cos(bet)-1)/2)*        exp(-1i*(+1)*gam), ...
    exp(-1i*(-2)*alp)*  (+sqrt(3/8)*sin(bet)^2)*           exp(-1i*( 0)*gam), ...
    exp(-1i*(-2)*alp)*  (+sin(bet)*(cos(bet)+1)/2)*        exp(-1i*(-1)*gam), ...
    exp(-1i*(-2)*alp)*  (+cos(bet/2)^4)*                   exp(-1i*(-2)*gam)];

end

% Consistency enforcement
function grumble(alp,bet,gam)
if (~isnumeric(alp))||(~isnumeric(bet))||(~isnumeric(gam))
    error('all inputs must be numeric.');
end
if (~isreal(alp))||(~isreal(bet))||(~isreal(gam))
    error('all inputs must be real.');
end
if (numel(alp)~=1)||(numel(bet)~=1)||(numel(gam)~=1)
    error('all inputs must have one element.');
end
end

% There is nothing of any importance in life - except how well you do
% your work. Nothing. Only that. Whatever else you are, will come from
% that. It's the only measure of human value. All the codes of ethics
% they'll try to ram down your throat are just so much paper money put
% out by swindlers to fleece people of their virtues. The code of com-
% petence is the only system of morality that's on a gold standard. 
% When you grow up, you'll know what I mean.
%
% Ayn Rand, "Atlas Shrugged"

