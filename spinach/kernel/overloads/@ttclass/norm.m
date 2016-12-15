% Computes the norm of the matrix represented by a tensor train. Syntax:
%
%                       ttnorm=norm(ttrain,norm_type)
%
%       norm_type=1         returns the 1-norm
%       norm_type=inf       returns the inf-norm
%       norm_type=2         returns the 2-norm
%       norm_type='fro'     returns the Frobenius norm
%
% Note: norms other than Frobenius norm are expensive for tensor trains.
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function ttnorm=norm(ttrain,norm_type)

% The default option is TYP=1
if (nargin == 1); norm_type=1; end

% Compute the norm
switch norm_type
    case 'fro'
        % Frobenius norm
        ttrain=ttclass_sum(ttrain);
        ttrain=ttclass_ort(ttrain,-1);
        ttnorm=abs(ttrain.coeff)*norm(ttrain.cores{1,1}(:));
    case 1
        % Maximum absolute column sum
        error('1-norm is not yet implemented for ttclass');
    case inf
        % Maximum absolute row sum
        error('inf-norm is not yet implemented for ttclass');
    case 2
        % Maximum absolute eigenvalue
        error('2-norm is not yet implemented for ttclass');
    otherwise
        % Complain and bomb out
        error('unrecognized norm type.');
end

end

% The most exciting phrase to hear in science, the one that heralds new
% discoveries, is not "eureka!" but rather "hmm....that's funny..."
%
% Isaac Asimov
%
%
% Contrary to what Asimov says, the most exciting phrase in science, the
% one that heralds new discoveries, is not "eureka!" or "that's funny...",
% it's "your research grant has been approved".
%
% John Alejandro King

