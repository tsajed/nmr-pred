% Adds corresponding fields of two structures. Nested structu-
% res are processed recursively. Syntax:
%
%          output_structure=plus(structure1,structure2)
%
% luke.edwards@ucl.ac.uk
% i.kuprov@soton.ac.uk

function str3=plus(str1,str2)

% Decide how to proceed
if isstruct(str1)&&isstruct(str2)
    
    % Get the field names
    fnames1=fieldnames(str1); fnames2=fieldnames(str2);
    
    % Check topology
    if (numel(fnames1)~=numel(fnames2))||...
       (~isempty(setdiff(fnames1,fnames2)))
        error('structure topology mismatch.');
    end
    
    % Loop over field names
    for n=1:length(fnames1)
        
        % Recursive call for each field name
        str3.(fnames1{n})=str1.(fnames1{n})+str2.(fnames1{n});
                      
    end
   
else
    
    % Complain and bomb out
    error('structure topology mismatch.');
    
end

end

% He was the sort of person who stood on mountaintops during
% thunderstorms in wet copper armour shouting "All the Gods 
% are bastards!"
%
% Terry Pratchett

