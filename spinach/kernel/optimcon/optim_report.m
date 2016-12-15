% Writes a log message to the console or an ACSII file. Syntax:
%
% <http://spindynamics.org/wiki/index.php?title=Optim_report.m>

function optim_report(output,report_string,depth,close_file)

% close all open files if selected
if nargin==4 && islogical(close_file) && close_file, 
    fclose('all'); 
    return; 
end

% Ignore the call if the system is hushed
if ~strcmp(output,'hush')
    
    % Validate the input
    grumble(output,report_string);
    
    % Compose the prefix
    call_stack=dbstack;
    for n=1:numel(call_stack)
        call_stack(n).name=[call_stack(n).name ' > '];
    end
    
    if ~exist('depth','var'), depth=0; 
    elseif depth>= numel(call_stack), depth=numel(call_stack)-1; end
    
    prefix_string=[call_stack(end:-1:2+depth).name];
    prefix_string=prefix_string(1:(end-3));
    
    % Fix empty prefixes
    if isempty(prefix_string)
        prefix_string=' ';
    end
    
    % Roll the prefix
    if numel(prefix_string)<50
        prefix_string=[prefix_string blanks(50)];
        prefix_string=prefix_string(1:50);
    else
        prefix_string=['...' prefix_string((end-46):end)];
    end
    
    % Add prefix to the report string
    report_string=['[' prefix_string ' ]  ' report_string];
    
    % Send the report string to the output
    fprintf(output,'%s\n',report_string);
    
end

end

% Consistency enforcement
function grumble(output,report_string)
if isempty(output)
    error('optim_system.sys.output field must exist.');
end
if ((~isa(output,'double'))&&(~isa(output,'char')))||...
   (isa(output,'char')&&(~strcmp(output,'hush')))
    error('optim_system.sys.output must be either ''hush'' or a file ID.');
end
if ~ischar(report_string)
    error('report_string must be a string.');
end
end
