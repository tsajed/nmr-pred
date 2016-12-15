% Number of nonzeros in all cores of a tensor train. Syntax: 
%
%                     answer=nnz(ttrain)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function answer=nnz(ttrain)

% Count the non-zeros
answer=sum(sum(cellfun(@nnz,ttrain.cores)));

end

% Men have always been and forever would remain silly victims of lies and
% self-deceit in politics, until they learn to see, behind any moral, re-
% ligious, political or social statements, proclamations and promises the
% interests of specific social classes.
%
% Vladimir Lenin

