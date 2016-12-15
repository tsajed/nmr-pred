% Multiplies all entries of a structure by a user-specified
% scalar. Nested structures are processed recursively. Syntax:
%
%      output_structure=mtimes(scalar,input_structure)
%
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk

function str_out=mtimes(scalar,str_in)

% Decide how to proceed
if isscalar(scalar)
    
    % Get the field names
    fnames=fieldnames(str_in);
    
    % Loop over field names
    for n=1:numel(fnames)
        
        % Recursive call for each field name
        str_out.(fnames{n})=scalar*str_in.(fnames{n});
                    
    end
    
else
    
    % Complain and bomb out
    error('the first argument must be a scalar.');
    
end

end

% Arthur Dent:  What happens if I press this button?
% Ford Prefect: I wouldn't --
% Arthur Dent:  Oh.
% Ford Prefect: What happened?
% Arthur Dent:  A sign lit up, saying "please do not press this button again".
%
% Douglas Adams

