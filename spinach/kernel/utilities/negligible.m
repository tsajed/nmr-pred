% An aux function determining whether a given object deserves attention
% given the tolerance specified. Used in the internal decision making 
% performed by Spinach kernel functions. Syntax:
%
%                answer=negligible(object,tolerance)
%
% Parameters:
%
%        object - the object whose significance is to be assessed
%
%     tolerance - significance threshold
%
% The function returns a logical value. If the object cannot be asses-
% ed based on the current heuristics, it is deemed non-negligible.
%
% i.kuprov@soton.ac.uk

function answer=negligible(object,tolerance)

answer=~significant(object,tolerance);

end

% It's too bad that stupidity isn't painful.
%
% Anton Szandor LaVey

