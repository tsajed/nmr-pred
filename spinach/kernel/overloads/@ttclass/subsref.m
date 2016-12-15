% Dot and bracket property specifications for the tensor train class.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function answer=subsref(ttrain,reference)

switch reference(1).type
    
    % Methods and properties
    case '.'
        
        % Return the output requested
        switch reference(1).subs 
            case 'ncores',    answer=size(ttrain.cores,1);
            case 'ntrains',   answer=size(ttrain.cores,2);
            case 'sizes',     answer=sizes(ttrain); 
            case 'ranks',     answer=ranks(ttrain);           
            case 'coeff',     answer=ttrain.coeff;
            case 'cores',     answer=ttrain.cores;
            case 'tolerance', answer=ttrain.tolerance;
            otherwise,        error(['unknown field reference ',reference(1).subs]);
        end
    
    % Core extraction
    case '{}'
        
        % Return the core requested
        if numel(reference(1).subs)~=2
            error('exactly two indices required to extract tt{j,k} kernel.');
        else
            answer=ttrain.cores{reference(1).subs{1},reference(1).subs{2}};
        end

    % Matrix element extraction
    case '()'
        
        % Start with zero
        answer=0;
        
        % Convert indices
        if numel(reference(1).subs)~=2
            error('exactly two indices required to evaluate tt(j,k) element');
        elseif (numel(reference(1).subs{1})==1)&&(numel(reference(1).subs{2})==1)
            siz=sizes(ttrain);
            ind=ttclass_ind2sub(siz(:,1),reference(1).subs{1});
            jnd=ttclass_ind2sub(siz(:,2),reference(1).subs{2});
        else
            error('advanced indexing is not implemented for tensor trains.');
        end
        
        % Multiply up the tensor train
        [ncores,ntrains]=size(ttrain.cores);
        for n=1:ntrains
            x=ttrain.cores{ncores,n}(:,ind(ncores),jnd(ncores),:);
            for k=(ncores-1):(-1):1
                x=ttrain.cores{k,n}(:,ind(k),jnd(k),:)*x;
            end
            answer=answer+ttrain.coeff(n)*x(1,1);
        end
        
    otherwise
        error('unknown subscript reference type.');
        
end

% Allow nested indexing
if numel(reference)>1
    answer=subsref(answer,reference(2:end));
end

end

% "A cactus is a very disappointed cucumber."
%
% A Russian saying

