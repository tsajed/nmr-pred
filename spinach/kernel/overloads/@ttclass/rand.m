% Generates tensor train matrix with random entries, same physical index
% topology at the tensor train supplied, and specified TT-ranks. Syntax:
%
%                         tt=rand(tt,ttrank)
%
% d.savostyanov@soton.ac.uk
% i.kuprov@soton.ac.uk

function tt=rand(tt,ttrank)

% Read tensor train sizes
[d,~]=size(tt.cores); sz=sizes(tt);

% Reallocate cores
tt.cores=cell(d,1);

% Fill the cores by random elements
tt.cores{1,1}=rand(1,sz(1,1),sz(1,2),ttrank);
for k=2:d-1
    tt.cores{k,1}=rand(ttrank,sz(k,1),sz(k,2),ttrank);
end
tt.cores{d,1}=rand(ttrank,sz(d,1),sz(d,2),1);

% Unit coefficient and zero tolerance
tt.coeff=1; tt.tolerance=0;

end

% The problem with inviting professors to dinner parties is that they're
% used to talking for a minimum of 50 minutes at a time.
%
% Anonymous philosophy professor

