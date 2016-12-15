% Returns the image encoded within the Fokker-Planck vector by the
% user-specified Liouville space coil state. Syntax:
%
%                    phan=fpl2phan(rho,coil,dims)
%
% Parameters:
%
%     rho    - state vector in Fokker-Planck space
%
%    coil    - observable state vector in Liouville space
%
%    dims    - spatial dimensions of the Fokker-Planck
%              problem
%
% i.kuprov@soton.ac.uk

function phan=fpl2phan(rho,coil,dims)

% Expose the spin dimension
rho=reshape(rho,[numel(coil) prod(dims)]);

% Compute the observable
phan=coil'*rho;

% Reshape as needed
phan=reshape(phan,dims);

end

% Malcolm Levitt's most fearsome battle cry, the one that sends the chill
% down the spines of any committee and makes marrow freeze in their bones
% in expectation of what is to come, is "Erm, excuse me..."

