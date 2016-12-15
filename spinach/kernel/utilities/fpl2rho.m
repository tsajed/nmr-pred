% Integrates over the spatial degrees of freedom and returns the
% average spin state vector across the sample. Syntax:
%
%                     rho=fpl2rho(rho,dims)
%
% Parameters:
%
%            rho   - Fokker-Planck state vector
%
%           dims   - spatial dimensions of the 
%                    Fokker-Planck problem
%
% i.kuprov@soton.ac.uk

function rho=fpl2rho(rho,dims)

% Expose the spin dimension
rho=reshape(rho,[numel(rho)/prod_dims prod(dims)]);

% Sum over the spatial coordinates
rho=sum(rho,2)/prod(dims);

end

% I would like to die on Mars, just not on impact.
%
% Elon Musk

