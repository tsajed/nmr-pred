% Chirp pulse waveform with a sine bell amplitude envelope in amplitude-
% frequency coordinates. Syntax:
%
%   [ampls,freqs]=chirp_pulse_af(npoints,duration,bandwidth,smfactor)
%
% with the following parameters:
%
%     npoints     - number of discretization points in
%                   the waveform
%
%    duration     - pulse duration, seconds
%
%    bandwidth    - chirp sweep bandwidth around 
%                   zero frequency, Hz
%
%    smfactor     - sine power to be used for the
%                   amplitude envelope
%
% i.kuprov@soton.ac.uk
% jeannicolas.dumez@cnrs.fr

function [ampls,freqs]=chirp_pulse_af(npoints,duration,bandwidth,smfactor)

% Compute amplitudes
ampls=1-abs(sin(linspace(-pi/2,pi/2,npoints)).^smfactor);

% Compute frequencies
freqs=linspace(-bandwidth/2,bandwidth/2,npoints);

% Calibrate amplitudes
ampls=2*pi*sqrt(bandwidth/duration)*ampls;

end

% My suggestion was quite simple: put that needed code number in a
% little capsule, and then implant that capsule right next to the
% heart of a volunteer. The volunteer would carry with him a big,
% heavy butcher knife as he accompanies the President. If ever the
% President wanted to fire nuclear weapons, the only way he could
% do so would be for him first, with his own hands, to kill one hu-
% man being. The President says, "George, I'm sorry but tens of mil-
% lions must die." He has to look at someone and realize what death
% is - what an innocent death is. Blood on the White House carpet. 
% Its reality brought home.
%
% When I suggested this to friends in the Pentagon they said, "My
% God, that's terrible. Having to kill someone would distort the
% President's judgment. He might never push the button."
%
% Roger Fisher

