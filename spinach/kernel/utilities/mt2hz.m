% Converts hyperfine couplings from milliTesla to Hz 
% (linear frequency). Syntax:
%
%                hfc_hz=mt2hz(hfc_mt)
%
% Arrays of any dimension are supported.
%
% i.kuprov@soton.ac.uk

function hfc_hz=mt2hz(hfc_mt)

if isnumeric(hfc_mt)&&isreal(hfc_mt)
    hfc_hz=1e7*2.802495365*hfc_mt;
else
    error('the input argument must be an array of real numbers.');
end
    
end

% You know that I write slowly. This is chiefly because I am never
% satisfied until I have said as much as possible in a few words, and
% writing briefly takes far more time than writing at length.
%
% Carl Friedrich Gauss

