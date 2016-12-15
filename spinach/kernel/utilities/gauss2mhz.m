% Converts hyperfine couplings from Gauss to MHz (linear
% frequency). Syntax:
%
%               hfc_mhz=gauss2mhz(hfc_gauss)
%
% Arrays of any dimensions are supported.
%
% i.kuprov@soton.ac.uk

function hfc_mhz=gauss2mhz(hfc_gauss)

if isnumeric(hfc_gauss)&&isreal(hfc_gauss)
    hfc_mhz=2.802495365*hfc_gauss;
else
    error('the argument must be an array of real numbers.');
end
    
end

% A celibate clergy is an especially good idea, because it tends to
% suppress any hereditary propensity toward fanaticism.
%
% Carl Sagan

