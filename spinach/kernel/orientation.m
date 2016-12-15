% Anisotropic part of the Hamiltonian for a specific spin system
% orientation. Syntax:
%
%                  H=orientation(Q,euler_angles)
% 
% Arguments:
%
%   Q             -  rotational basis as returned by hamiltonian
%                    function.
%
%   euler_angles  -  a 1x3 vector or a vertical stack of 1x3 
%                    vectors specifying Euler angles (radians)
%                    relative to the input orientation.
%
% Output:
%
%     H    -    anisotropic part of the Hamiltonian commutation
%               superoperator (or a cell array thereof) for the 
%               specified Euler angles.
%
% This function may be used in both Hilbert and Liouville space
% because the H -> [H, ] adjoint map is linear.
%
% i.kuprov@soton.ac.uk

function H=orientation(Q,euler_angles)

% Check consistency
grumble(Q,euler_angles);

% Determine problem dimension
n_orientations=size(euler_angles,1); 

% Preallocate the result array
H=cell(n_orientations,1);
for n=1:n_orientations
    H{n}=0*Q{1,1};
end

% Compute Wigner matrices
D=cell(n_orientations,1);
for n=1:n_orientations
    D{n}=euler2wigner(euler_angles(n,1),...
                      euler_angles(n,2),...
                      euler_angles(n,3));
end

% Compute Hamiltonian operators
for n=1:n_orientations
    for k=1:5
        for m=1:5
            H{n}=H{n}+(D{n}(k,m))*Q{k,m};
        end
    end
    H{n}=(H{n}+H{n}')/2;
end

% If there is only one operator to return, remove the cell
if numel(H)==1, H=H{1}; end

end

% Consistency enforcement
function grumble(Q,euler_angles)
if (~iscell(Q))||any(size(Q)~=[5 5])
    error('Q parameter must be a 5x5 cell array of matrices.');
end
if (~isnumeric(euler_angles))||(size(euler_angles,2)~=3)
    error('Euler angles must be supplied as a row vector or a vertical stack thereof.')
end
end

% 1f y0u c4n r34d 7h15, y0u r34||y n33d 70 637 |41d

