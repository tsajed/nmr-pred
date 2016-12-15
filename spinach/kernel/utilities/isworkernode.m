% Returns true if executed inside a parfor or spmd block.
%
% i.kuprov@soton.ac.uk

function answer=isworkernode()

    % Best way I could find
    answer=~isempty(getCurrentTask());
    
end

% In the beginning the Universe was created. This has 
% made a lot of people very angry and been widely re-
% garded as a bad move.
%
% Douglas Adams, "Hitchhiker's Guide to the Galaxy"

