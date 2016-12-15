% Chirp pulse waveform with a sine bell amplitude envelope in Cartesian
% coordinates. Syntax:
%
%       [Cx,Cy]=chirp_pulse_xy(npoints,duration,bandwidth,smfactor)
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

function [Cx,Cy]=chirp_pulse_xy(npoints,duration,bandwidth,smfactor)

% Compute timing grid
time_grid=linspace(-0.5,0.5,npoints);

% Compute amplitudes
amplitudes=1-abs(sin(pi*time_grid).^smfactor);

% Compute phases
phases=pi*duration*bandwidth*(time_grid.^2);

% Calibrate amplitudes
amplitudes=2*pi*sqrt(bandwidth/duration)*amplitudes;

% Transform into Cartesian representation
[Cx,Cy]=polar2cartesian(amplitudes,phases);

end

% "No one realized that the pumps that delivered fuel to the 
% emergency generators were electric."
%
% Angel Feliciano, representative of Verizon
% explaining why Verizon's backup power failed
% during the August 13, 2003 blackout, causing
% disruption to the 911 service.

