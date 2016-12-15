% Number of elements in the matrix represented by a tensor train.
% Syntax:
%                        n=ttclass_numel(ttrain)
%
% For large spin systems it may be too large to fit the integer.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function n=ttclass_numel(ttrain)

% Compute the number of elements
if isa(ttrain,'ttclass')
    n=prod(prod(sizes(ttrain)));
else
    error('this function only applies to tensor trains.');
end

% Check for infinities
if n>intmax
    error('the number of elements exceeds Matlab''s intmax.');
end

end

% If it had been possible to build the tower of Babel without
% ascending it, the work would have been permitted. 
%
% Franz Kafka

