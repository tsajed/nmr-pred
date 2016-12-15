% Wipes the diagonal using the brush with the specified dimensions.
%
% i.kuprov@soton.ac.uk

function spec=killdiag(spec,brush_dim)

% Loop over the column index
for n=1:size(spec,2)
    
    % Find the row index
    k=n*size(spec,1)/size(spec,2);
    
    % Find the row index extents
    k=floor(k-brush_dim/2):ceil(k+brush_dim/2);
    
    % Avoid array boundaries
    k(k<1)=[]; k(k>size(spec,1))=[];
    
    % Zero the elements
    spec(k,n)=0;
    
end

end

% "I do wish we could chat longer, but... I'm having
%  an old friend for dinner."
%
% Hannibal Lecter

