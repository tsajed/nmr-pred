% Converts Hartree energy units into J/mol. Syntax:
%
%           energy=hartree2joule(energy)
%
% Arrays of any dimension are supported.
%
% i.kuprov@soton.ac.uk

function energy=hartree2joule(energy)

% Perform the conversion
energy=2625499.62*energy;

end

% "We only need to be lucky once. You need to be
%  lucky every time."
%
% The IRA to Margaret Thatcher, after 
% a failed assassination attempt.

