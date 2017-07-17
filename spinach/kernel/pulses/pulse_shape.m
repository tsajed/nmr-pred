% Shaped pulse waveforms. Add your own if necessary. Syntax:
%
%             waveform=pulse_shape(pulse_name,npoints)
%
% Parameters:
%
%     pulse_name - the name of the pulse
%
%     npoints    - number of points in the pulse
%
% i.kuprov@soton.ac.uk

function waveform=pulse_shape(pulse_name,npoints)

% Check consistency
grumble(pulse_name,npoints);

% Choose the shape
switch pulse_name
    
    case 'gaussian'
        
        time_grid=linspace(-2,2,npoints);
        waveform=normpdf(time_grid);
        
    case 'rectangular'
        
        waveform=ones(1,npoints);
        
    otherwise
        
        % Complain and bomb out
        error('unknown pulse name.');
        
end

end

% Consistency enforcement
function grumble(pulse_name,npoints)
if ~ischar(pulse_name)
    error('pulse_name parameter must be a character string.');
end
if (numel(npoints)~=1)||(~isnumeric(npoints))||(~isreal(npoints))||...
   (npoints<1)||(mod(npoints,1)~=0)
    error('npoints parameter must be a positive real integer greater than 1.');
end
end

% Let me start with a parable. It concerns an Eastern European country
% whose parliament was considering a total smoking ban. In response, a
% consortium of tobacco companies demonstrated that the savings made in
% healthcare as a result of the decline in smoking-related diseases were
% chicken feed besides the reduced payout in pensions as the result of
% premature death - not to mention the fiscal increment from the habit.
%
% Stewart Dakers
