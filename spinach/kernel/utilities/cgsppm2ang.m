% Converts magnetic susceptibility from the cgs-ppm (aka cm^3/mol) units
% quoted by quantum chemistry packages into Angstrom^3 units required by
% Spinach pseudocontact shift functionality. Syntax:
%
%                        ang=cgsppm2ang(cgsppm)
%
% Arrays of any dimension are supported.
%
% i.kuprov@soton.ac.uk

function ang=cgsppm2ang(cgsppm)

% Check consistency
grumble(cgsppm);

% Do the calculation
ang=4*pi*1e18*cgsppm/6.02214129e23;

end

% Consistency enforcement
function grumble(cgsppm)
if ~isnumeric(cgsppm)
    error('input must be numeric.');
end
end

% "What is this thing, anyway?" said the Dean, inspecting the implement in
% his hands. "It's called a shovel," said the Senior Wrangler. "I've seen
% the gardeners use them. You stick the sharp end in the ground. Then it
% gets a bit technical."
%
% Terry Pratchett

