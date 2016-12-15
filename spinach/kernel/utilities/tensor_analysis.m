% Returns diagnostic information about an interaction tensor. Syntax:
%
%        [eigvals,dcm,iso]=tensor_analysis(spin_system,tensor)
%
% Parameters:
%
%           tensor - a 3x3 Cartesian interaction tensor matrix 
%
% Outputs:
%
%          eigvals - eigenvalues of the tensor
%
%              dcm - directional cosine matrix for the orientation
%                    of the tensor eigenframe relative to the fra-
%                    me it's been supplied in.
%
%              iso - isotropic part of the tensor
%
% Note: directional cosine matrix is not unique. Although the function
%       makes some effort at consistency, it cannot be guaranteed.
%
% i.kuprov@soton.ac.uk

function [eigvals,dcm,iso]=tensor_analysis(spin_system,tensor)

% Check consistency
grumble(tensor);

% Get the eigensystem
[dcm,eigvals]=eig(symmetrize(spin_system,full(tensor)));

% Sort the eigensystem
[~,index]=sort(abs(diag(eigvals)),1,'ascend');

% Rearrange the eigenvalues
eigvals=diag(eigvals);
eigvals=eigvals(index);

% Rearrange the eigenvectors
dcm=dcm(:,index);

% Kill the inversion component
dcm=dcm*det(dcm);

% Prefer upper half-space for positive directions
if dcm(3,3)<0
    dcm(:,3)=-dcm(:,3);
    dcm(:,1)=-dcm(:,1);
end

% Compute the isotropic part
iso=mean(eigvals);

end

% Consistency enforcement
function grumble(tensor)
if (~isnumeric(tensor))||(~ismatrix(tensor))||any(size(tensor)~=[3 3])||(~isreal(tensor))
    error('interaction tensor must be a real 3x3 matrix.');
end
end

% He who dares not offend cannot be honest.
%
% Thomas Paine

