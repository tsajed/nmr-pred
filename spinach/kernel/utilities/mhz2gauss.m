% Converts hyperfine couplings from MHz (linear frequency)
% to Gauss. Syntax:
%
%               hfc_gauss=mhz2gauss(hfc_mhz)
%
% Arrays of any dimensions are supported.
%
% i.kuprov@soton.ac.uk

function hfc_gauss=mhz2gauss(hfc_mhz)

if isnumeric(hfc_mhz)&&isreal(hfc_mhz)
    hfc_gauss=hfc_mhz/2.802495365;
else
    error('the argument must be an array of real numbers.');
end
    
end

% "Trillian had come to suspect that the main reason [Zaphood] had 
% had such a wild and successful life was that he never really un-
% derstood the significance of anything he did."
%
% Douglas Adams, "The Hitchhiker's Guide to the Galaxy"

