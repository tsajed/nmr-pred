% Generates spherical tensor expansions for single-spin Zeeman
% basis operators in Hilbert space. Syntax:
%
%                stexp=zeeman2sphten(matrix,type)
%
% where matrix is a single-spin matrix written in the Zeeman ba-
% sis in Hilbert space, type denote the matrix is a density mat-
% rix or a operator, and stexp is a cell array with the first 
% element of each row giving the operator name and the second 
% element being the corresponding expansion coefficient. The se-
% cond argument should be set to 'state' to get a state vector
% expansion and to 'oper' to get an operator expansion.
%
% lipingyang87@gmail.com
% i.kuprov@soton.ac.uk

function stexp=zeeman2sphten(matrix,type)

% Validate the input
grumble(type);

% Get the IST operators
I=irr_sph_ten(numel(diag(matrix)));

% Preallocate the output
stexp=cell(numel(I),2);

% Run the expansion
for n=1:numel(I)
    
    % Build operator name
    [L,M]=lin2lm(n-1);
    stexp{n,1}=['T' num2str(L) ',' num2str(M)];
    
    % Choose correct normalisation
    switch type
        
        case 'state'
        
            % Get the coefficient
            stexp{n,2}=trace(I{n}'*matrix)/sqrt(trace(I{n}'*I{n}));
        
        case 'oper'
            
            % Get the coefficient
            stexp{n,2}=trace(I{n}'*matrix)/trace(I{n}'*I{n});            
        
        otherwise
            
            % Complain and bomb out
            error('The input type is wrong.');
            
    end
    
end

end

% Input validation function
function grumble(type)
if ~ischar(type)
    error('the variable type must be a character string.');
end
if ~ismember(type,{'state','oper'})
    error('valid values for type are ''state'' and ''oper''.');
end
end

