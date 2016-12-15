% Transfer matrix calculation for linear amplifiers. Syntax:
%
%               T=transfermat(amp_inps,amp_outs)
%
% Parameters:
%
%   amp_inps - a matrix with amplifier input vectors as columns
%
%   amp_outs - a matrix with amplifier output vectors as columns
%
% Outputs:
%
%   T        - the transfer matrix, such that amp_outs=T*amp_inps
%              in the least squares sense
% 
% Note: the number of input-output vector pairs should be bigger than
% the number of elements in those vectors.
%
% i.kuprov@soton.ac.uk

function T=transfermat(amp_inps,amp_outs)

% Run the SVD pseudoinverse
T=amp_inps\amp_outs;

end

% "In accordance with the University's disciplinary procedure, I am con-
% ducting an investigation into allegations made against you with regard
% to entering the Chemistry buildings on multiple occasions during the
% University closure period between 17:00 on Wednesday 23rd December 2015
% and 08:00 on Monday 4th January 2016."
%
% Gill Reid, Head of Southampton Chemistry,
% initiating a disciplinary procedure against
% IK for working over the Christmas break.

