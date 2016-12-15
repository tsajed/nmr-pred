% Normalized Gaussian function in magnetic resonance notation. Syntax:
%
%                          g=gaussfun(x,fwhm)
%
% where fwhm is the full width at half-maximum.
%
% i.kuprov@soton.ac.uk

function y=gaussfun(x,fwhm)

% Compute standard deviation
sigma=fwhm/(2*sqrt(2*log(2)));

% Compute the Gaussian
y=(1/(sigma*sqrt(2*pi)))*exp(-(x.^2)/(2*sigma^2));

end

% Fifty years ago the back streets of Leningrad
% have taught me one lesson: when a fight is un-
% avoidable, punch first.
%
% Vladimir Putin

