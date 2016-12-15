% Reads wavegauss.pk files. Arguments:
%
%   filename   -   a string containing the name of the file
%
%   npoints    -   waveform upsampling or downsampling is
%                  performed to this number of points
%
% Note: phases are returned in degrees, amplitudes in percent.
% 
% i.kuprov@soton.ac.uk

function [amplitudes,phases]=read_wave(filename,npoints)

% Read the file
wavefile=fopen(filename,'r');
waveform=textscan(wavefile,'%f, %f','CommentStyle','##');
fclose(wavefile);

% Build complex waveform
waveform=waveform{1}.*exp(1i*pi*waveform{2}/180);

% Resample the waveform
waveform=interp1(linspace(0,1,numel(waveform)),real(waveform),linspace(0,1,npoints),'pchip')+...
      1i*interp1(linspace(0,1,numel(waveform)),imag(waveform),linspace(0,1,npoints),'pchip');

% Convert back into amplitudes and phases
amplitudes=abs(waveform);
phases=90-(180/pi)*atan2(real(waveform),imag(waveform));

end

% I condemn Christianity... It is, to me, the greatest of all imaginable
% corruptions. [...] To breed in humans a self-contradiction, an art of
% self-pollution, a will to lie at any price, an aversion and contempt for
% all good and honest instincts! [...] The beyond as the will to deny all
% reality; the cross as the distinguishing mark of the most subterranean
% conspiracy ever heard of -- against health, beauty, well-being, intel-
% lect... against life itself.
%
% I call Christianity the one great curse, the one great intrinsic depra-
% vity, the one great instinct of revenge, for which no means are veno-
% mous enough, or secret, subterranean and small enough - I call it the 
% one immortal blemish upon the human race...
%
% Friedrich Nietzsche

