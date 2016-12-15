% Approximates a given tensor train with lower TT-ranks. Syntax:
%
%                         ttrain=shrink(ttrain)
%
% Returns a single tensor train with right-to-left orthogonalisation.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function ttrain=shrink(ttrain)

% Read train sizes
[d,~]=size(ttrain.cores);

% Summation
ttrain=ttclass_sum(ttrain);

% Right-to-left orthogonalization
ttrain=ttclass_ort(ttrain);

% Check the norm and escape if the object is zero
nrm=ttrain.coeff*norm(ttrain.cores{d,1}(:));
if nrm==0, ttrain=0*unit_like(ttrain); return; end

% Truncation
ttrain=ttclass_svd(ttrain);

% Convert to a scalar if appropriate
if all(all(cellfun(@(x)size(x,2),ttrain.cores)==1))&&...
   all(all(cellfun(@(x)size(x,3),ttrain.cores)==1))
    ttrain=full(ttrain);
end
    
end

% "It's a tough life, being small and delicious."
%
% A Russian saying

