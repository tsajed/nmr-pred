% Writes a log message to the console or an ACSII file. Syntax:
%
%                     report(spin_system,report_string)
%
% where report_string is a character string. A newline symbol at the end
% of the string is not necessary - it is added by the function.
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function report(spin_system,report_string)

% Catch single-argument calls
if nargin==1
    error('console reporting function requires two arguments.');
end

% Ignore the call if the system is hushed
if ~strcmp(spin_system.sys.output,'hush')
    
    % Validate the input
    grumble(spin_system,report_string);
    
    % Compose the prefix
    call_stack=dbstack;
    for n=1:numel(call_stack)
        call_stack(n).name=[call_stack(n).name ' > '];
    end
    prefix_string=[call_stack(end:-1:2).name];
    prefix_string=prefix_string(1:(end-3));
    
    % Fix empty prefixes
    if isempty(prefix_string)
        prefix_string=' ';
    end
    
    % Roll the prefix
    if numel(prefix_string)<50
        prefix_string=pad(prefix_string,50);
    else
        prefix_string=['...' prefix_string((end-46):end)];
    end
    
    % Add prefix to the report string
    report_string=['[' prefix_string ' ]  ' report_string];
    
    % Send the report string to the output, ignoring impossible writes
    try fprintf(spin_system.sys.output,'%s\n',report_string); end %#ok<TRYNC>
    
end

end

% Consistency enforcement
function grumble(spin_system,report_string)
if (~isfield(spin_system,'sys'))||~isfield(spin_system.sys,'output')
    error('spin_system.sys.output field must exist.');
end
if ((~isa(spin_system.sys.output,'double'))&&(~isa(spin_system.sys.output,'char')))||...
   (isa(spin_system.sys.output,'char')&&(~strcmp(spin_system.sys.output,'hush')))
    error('spin_system.sys.output must be either ''hush'' or a file ID.');
end
if ~ischar(report_string)
    error('report_string must be a string.');
end
end

% All parts should go together without forcing. You must remember that the
% parts you are reassembling were disassembled by you. Therefore, if you
% can't get them together again, there must be a reason. By all means, do
% not use a hammer.
%
% IBM Manual, 1925 

