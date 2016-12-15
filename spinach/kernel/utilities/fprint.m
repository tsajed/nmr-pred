% 2D spectrum integral fingerprinting utility. Breaks the spectrum
% down into bins of user-specified size and returns the integrals
% over those bins. Syntax:
%
%                     fp=fprint(spectrum,nbins)
%
% Parameters:
%
%      spectrum    - the matrix of a real or complex 2D spectrum
%
%         nbins    - a vector with two integers, specifying the
%                    number of bins in each spectral dimension
%
% Note: spectrum matrix dimensions must be divisible by the corres-
% ponding element of nbins vector.
%
% i.kuprov@soton.ac.uk

function fp=fprint(spectrum,nbins)

% Decide problem dimensionality
switch ndims(spectrum)
    
    case 2
        
        % Preallocate the answer
        fp=zeros(nbins);
        
        % Determine bin boundaries
        bins_x=linspace(0,size(spectrum,1),nbins(1)+1);
        bins_y=linspace(0,size(spectrum,2),nbins(2)+1);
        
        % Loop over the bins
        for n=1:nbins(1)
            for k=1:nbins(2)
                
                % Compute bin integrals
                fp(n,k)=trapz(trapz(spectrum((bins_x(n)+1):bins_x(n+1),(bins_y(k)+1):bins_y(k+1))));
                
            end
        end
                
    otherwise
        
        error('not implemented');

end

% Never underestimate the bandwidth of a station wagon full of
% tapes hurtling down the highway.
%
% Andrew Tanenbaum

