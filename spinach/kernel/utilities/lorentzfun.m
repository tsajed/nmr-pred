% Normalized Lorentzian function in magnetic resonance notation. Syntax:
%
%                          g=lorentzfun(x,fwhm)
%
% where fwhm is the full width at half-maximum.
%
% i.kuprov@soton.ac.uk

function y=lorentzfun(x,fwhm)

% Get the width parameter
gamma=fwhm/2;

% Compute the Lorentzian
y=(1/(pi*gamma))*(1./(1+(x/gamma).^2));

end

% [...] the sexual revolution was the mockery of the century,
% because now women were giving to men for free what they used
% to have to marry us for.
%
% A feminist web site

