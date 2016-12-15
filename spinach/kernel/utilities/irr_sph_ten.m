% Returns a cell array of single-spin irreducible spherical tensor opera-
% tors T(k,m). A two-argument call
%
%                          T=irr_sph_ten(mult,k)
%
% where 'mult' is the multiplicity of the spin in question and 'k' is the
% irreducible spherical tensor rank required, returns a cell array of ten-
% sors of that rank in the order of decreasing projection. A single argu-
% ment call
%                           T=irr_sph_ten(mult)
%
% produces tensors of all ranks and concatenates them into a cell array in
% the order of ascending rank.
%
% The resulting spherical tensors are normalized in such a way as to obey
% the following commutation relation:
%
%                            [Lz,T_lm]=m*T_lm
%
% Note: operator normalization in spin dynamics is an old and thorny ques-
% tion. Many different conventions exist, but the only way to make the re-
% sulting formalism independent of the total spin quantum number is to im-
% pose identical commutation relations rather than equal matrix norms.
%
% i.kuprov@soton.ac.uk
% hannah.hogben@chem.ox.ac.uk

function T=irr_sph_ten(mult,k)

% Adapt to the input style
switch nargin
    
    case 1
        
        % Generate tensors of all ranks and put them into a cell array
        if isnumeric(mult)&&(numel(mult)==1)&&(mult>0)&&(mod(mult,1)==0)
            T=cell(mult^2,1);
            for n=0:(mult-1)
                T((n^2+1):((n+1)^2))=irr_sph_ten(mult,n);
            end
        else
            error('mult must be a non-negative integer.');
        end
        
    case 2
        
        % Generate spherical tensor operators of the specified rank
        if isnumeric(mult)&&(numel(mult)==1)&&(mod(mult,1)==0)&&...
           isnumeric(k)&&(numel(k)==1)&&(k>0)&&(k<mult)&&(mod(k,1)==0)
            
            % Get Pauli matrices
            L=pauli(mult);
            
            % Preallocate the cell array
            T=cell(2*k+1,1);
            
            % Get the top state
            T{1}=((-1)^k)*(2^(-k/2))*L.p^k;
            
            % Apply sequential lowering using Racah's commutation rule
            for n=2:(2*k+1)
                q=k-n+2; T{n}=(1/sqrt((k+q)*(k-q+1)))*(L.m*T{n-1}-T{n-1}*L.m);
            end
                
        elseif isnumeric(mult)&&(numel(mult)==1)&&(mod(mult,1)==0)&&...
               isnumeric(k)&&(numel(k)==1)&&(k==0)
            
            % For zero rank, return a unit matrix
            T={speye(mult)};
            
        else
            
            % For non-physical ranks, complain and bomb out
            error('mult must be a non-negative integer and k must be an integer from [0, mult-1] interval.');
            
        end
        
end

end

% I swear by my life and my love of it that I will never live for the sake
% of another man, nor ask another man to live for mine.
%
% Ayn Rand, "Atlas Shrugged"

