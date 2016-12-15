% Symmetrizes interaction tensors. This is required to avoid numerical
% instabilities associated with small asymmetries during the diagonali-
% zation process. Syntax:
%
%                      A=symmetrize(spin_system,A)
%
% i.kuprov@soton.ac.uk

function A=symmetrize(spin_system,A)

% Check consistency
grumble(A);

% If the asymmetry is significant, make a fuss
if norm(A-transpose(A))>spin_system.tols.inter_sym
    report(spin_system,'WARNING - significant asymmetry detected in a coupling tensor.');
    report(spin_system,['WARNING - symmetric part norm: ' num2str(norm((A+A')/2)/(2*pi)) ' Hz']);
    report(spin_system,['WARNING - antisymmetric part norm: ' num2str(norm((A-A')/2)/(2*pi)) ' Hz']);
    report(spin_system,'WARNING - the tensor has been symmetrized.');
end

% Symmetrize the tensor
A=(A+transpose(A))/2;

end

% Consistency enforcement
function grumble(tensor)
if (~isnumeric(tensor))||(~ismatrix(tensor))||any(size(tensor)~=[3 3])||(~isreal(tensor))
    error('interaction tensor must be a real 3x3 matrix.');
end
end

% "Smoking - NO HYDROGEN!"
%
% A safety warning on Anatole Abragam's door

