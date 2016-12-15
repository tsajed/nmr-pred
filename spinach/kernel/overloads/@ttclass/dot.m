% Dot product of TT representations of matrices. Syntax:
%
%                        c=dot(a,b)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function c=dot(a,b)

% Check consistency
grumble(a,b);

% Compute the product
c=ctranspose(a)*b;

end

% Consistency enforcement
function grumble(a,b)
if (~isa(a,'ttclass'))||(~isa(b,'ttclass'))
    error('both inputs should be tensor trains.')
end
if (size(a.cores,1)~=size(b.cores,1))||(~all(all(sizes(a)==sizes(b))))
    error('tensor train structures must be consistent.')
end
end

% Any product that needs a manual to work is broken.
%
% Elon Musk

