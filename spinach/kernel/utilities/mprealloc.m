% Preallocates an operator in the current basis. Syntax:
%
%             A=mprealloc(spin_system,nnzpc)
%
% Parameters:
%
%    nnzpc  -  expected number of non-zeros per column
%
% i.kuprov@soton.ac.uk

function A=mprealloc(spin_system,nnzpc)

% Check the input
if (~isnumeric(nnzpc))||(~isreal(nnzpc))||(~isscalar(nnzpc))||(mod(nnzpc,1)~=0)
    error('nnzpc parameter must be a positive real integer.');
end

% Do the math
switch spin_system.bas.formalism
    
    case 'sphten-liouv'
        
        % Create a zero Liouville space matrix operator
        problem_dim=size(spin_system.bas.basis,1);
        A=spalloc(problem_dim,problem_dim,nnzpc*problem_dim);
        
    case 'zeeman-hilb'
        
        % Create a zero Hilbert space matrix operator
        problem_dim=prod(spin_system.comp.mults);
        A=spalloc(problem_dim,problem_dim,nnzpc*problem_dim);
        
    case 'zeeman-liouv'
        
        % Create a zero Liouville space matrix operator
        problem_dim=prod(spin_system.comp.mults.^2);
        A=spalloc(problem_dim,problem_dim,nnzpc*problem_dim);
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% In 2006, Oxford's Magdalen College (where Erwin Schrodinger was a Fellow
% between 1933 and 1936) received a sum of money from a benefactor towards
% "increasing the art content of the College". A number of works were pre-
% sented for competition, among them a beautiful stone obelisk, called "Mo-
% nument to Knowledge" with the Schrodinger equation inscribed on it. The
% obelisk was rejected - the inscriber had missed the bar off Planck's con-
% stant. As the College Governing Body put it, the equation, as written,
% "would have exploded the stone it was inscribed upon".

