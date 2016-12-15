% Merges timing tables of pulse sequence channels. Both timing tables
% should be given as row arrays of the following structure:
%
%           A={A_1 A_2 ... A_n};   dtA=[dt_1 dt_2 ... dt_n];
%
% where A_k is a spin operator and dt_k is the time for which it acts.
% The sequence is assumed to be executed chronologically from left to
% right. The two event sequences can be aligned left, right or centre
% with respect to one another. Syntax:
%
%                  [C,dtC]=splice(A,dtA,B,dtB,alignment)
%
% i.kuprov@soton.ac.uk

function [C,dtC]=splice(A,dtA,B,dtB,alignment)

% Check consistency
grumble(A,dtA,B,dtB,alignment);

% Determine problem dimension
dim=size(A{1},1);

% Set the alignment type
switch alignment
    
    case 'left'
        
        % Find event edges
        tA_edges=[0 cumsum(dtA)]; tB_edges=[0 cumsum(dtB)];
        
        % Pad the shorter sequence with zero operators
        if max(tA_edges)<max(tB_edges), A=[A {spalloc(dim,dim,0)}]; end
        if max(tA_edges)>max(tB_edges), B=[B {spalloc(dim,dim,0)}]; end
        
        % Build new timing diagram
        tC_edges=unique([tA_edges tB_edges]);
        
        % Build new operator list
        C=cell(1,numel(tC_edges)-1);
        for n=1:numel(C)
            
            % Find values of A and B in the current period
            current_time=(tC_edges(n)+tC_edges(n+1))/2;
            A_part=A(find(tA_edges<current_time,1,'last'));
            B_part=B(find(tB_edges<current_time,1,'last'));

            % Assign the result
            if isempty(A_part), A_part={spalloc(dim,dim,0)}; end
            if isempty(B_part), B_part={spalloc(dim,dim,0)}; end
            C{n}=A_part{1}+B_part{1};
            
        end
        
        % Compute new time increments
        dtC=diff(tC_edges);
        
    case 'right'
        
        % Reverse the inputs
        A=flip(A); dtA=flip(dtA);
        B=flip(B); dtB=flip(dtB);
        
        % Perform left-aligned splice
        [C,dtC]=splice(A,dtA,B,dtB,'left');
        
        % Reverse the output
        C=flip(C); dtC=flip(dtC);
        
    case 'centre'
        
        % Find event edges
        tA_right=cumsum(dtA); tB_right=cumsum(dtB);
        
        % Pad the shorter sequence with zero operators
        if tA_right(end)<tB_right(end)
            dtA=[(tB_right(end)-tA_right(end))/2 dtA]; A=[{spalloc(dim,dim,0)} A];
        elseif tA_right(end)>tB_right(end)
            dtB=[(tA_right(end)-tB_right(end))/2 dtB]; B=[{spalloc(dim,dim,0)} B];
        end
        
        % Perform left-aligned splice
        [C,dtC]=splice(A,dtA,B,dtB,'left');
          
    otherwise
        
        % Complain and bomb out
        error('unknown alignment type.');
        
end

end

% Consistency enforcement
function grumble(A,dtA,B,dtB,alignment)
if (~iscell(A))||(~iscell(B))||isempty(A)||isempty(B)||(~isrow(A))||(~isrow(B))
    error('A and B arguments must be row cell arrays of matrices.');
end
if (~all(cellfun(@issparse,A)))||(~all(cellfun(@issparse,B)))
    error('all elements of A and B must be sparse.');
end
if (~isnumeric(dtA))||(~isnumeric(dtB))||(~isrow(dtA))||(~isrow(dtB))||any(dtA<=0)||any(dtB<=0)
    error('dtA and dtB arguments must be row vectors of positive numbers.');
end
if ~ischar(alignment)
    error('the alignment argument must be a character string.');
end
end

% The JEOL Student Prize, awarded at the annual conferences of the Royal Society
% of Chemistry ESR Group, had one winner in its history who got the Prize before
% she even started talking. As the laptop was being connected to the projection
% system and the usual technology gremlins refused to work, Petra Luders uncere-
% moniously pushed aside the lost-looking IT guy and got the thing to work in a
% few keystrokes on what appeared to be a Linux system. "And here's our winner"
% thought every Committee member with a quiet smile. The science and the presen-
% tation easily outshined every other competitor - Ms Luders got her JEOL Medal.

