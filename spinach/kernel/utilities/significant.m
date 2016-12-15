% An aux function determining whether a given object deserves attention
% given the tolerance specified. Used in the internal decision making 
% performed by Spinach kernel functions. Syntax:
%
%                answer=significant(object,tolerance)
%
% Parameters:
%
%        object - the object whose significance is to be assessed
%
%     tolerance - significance threshold
%
% The function returns a logical value. If the object cannot be asses-
% ed based on the current heuristics, it is deemed significant.
%
% i.kuprov@soton.ac.uk

function answer=significant(object,tolerance)

% Check the input
if (~isnumeric(tolerance))||(~isscalar(tolerance))||(~isreal(tolerance))||(tolerance<0)
    error('tolerance parameter must be a non-negative real number.');
end

% Decide on the significance
if ~isnumeric(object)
    
    % Non-numeric objects are significant
    answer=true();

elseif isempty(object)
    
    % Empty arrays are not significant
    answer=false();
    
elseif norm(object,1)<tolerance
    
    % Arrays with a small norm are not significant
    answer=false();

else
    
    % Everything else is significant
    answer=true();
    
end
       
end

% There is a beast in man that needs to be exercised, not exorcised.
%
% Anton Szandor LaVey

