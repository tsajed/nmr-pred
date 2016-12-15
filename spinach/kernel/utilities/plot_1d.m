% 1D plotting utility. Uses the same parameters structure as 1D pulse
% sequences. Input syntax:
%
%              plot_1d(spin_system,spectrum,parameters)
%
% Input parameters:
%
%    parameters.sweep              sweep width, Hz
%
%    parameters.spins              spin species, e.g. {'1H'}
%
%    parameters.offset             transmitter offset, Hz
%
%    parameters.axis_units         axis units ('ppm','Gauss',
%                                  'mT','T','Hz','kHz','MHz')
%
%    parameters.derivative         if set to 1, the spectrum is
%                                  differentiated before plotting
%
%    parameters.invert_axis        if set to 1, the frequency axis 
%                                  is inverted before plotting
%
% Any extra arguments given to this function will be passed to the built-in
% Matlab function plot().
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function plot_1d(spin_system,spectrum,parameters,varargin)

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spectrum,parameters);

% If a complex spectrum is received, plot both components
if ~isreal(spectrum)
    
    % Recursively plot the real component
    subplot(2,1,1); plot_1d(spin_system,real(spectrum),parameters,varargin{:});
    title('Real part of the complex spectrum.');
    
    % Recursively plot the imaginary component
    subplot(2,1,2); plot_1d(spin_system,imag(spectrum),parameters,varargin{:});
    title('Imaginary part of the complex spectrum.'); return;
    
end

% Get the axis
[ax,ax_label]=axis_1d(spin_system,parameters);

% Compute the derivative if necessary
if isfield(parameters,'derivative')&&parameters.derivative
    spectrum=fdvec(spectrum,5,1);
end

% Plot the spectrum
plot(ax,spectrum,varargin{:}); axis tight;

% Label the axis
xlabel(ax_label);

% Invert the axis if necessary
if isfield(parameters,'invert_axis')&&parameters.invert_axis
    set(gca,'XDir','reverse');
end

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if (~isfield(parameters,'offset'))&&(numel(parameters.sweep)==1)
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'axis_units')
    report(spin_system,'parameters.axis_units field not set, assuming ppm.');
    parameters.axis_units='ppm';
end
if ~isfield(parameters,'invert_axis')
    report(spin_system,'parameters.invert_axis field not set, assuming NMR tradition.');
    parameters.invert_axis=1;
end
end

% Consistency enforcement
function grumble(spectrum,parameters)
if (~isnumeric(spectrum))||(~isvector(spectrum))
    error('spectrum should be a vector of numbers.');
end
if (~isfield(parameters,'offset'))&&(numel(parameters.sweep)==1)
    error('offset should be specified in parameters.offset variable.');
end
if (numel(parameters.sweep)==2)&&isfield(parameters,'offset')
    error('offset should not be specified when spectral extents are given in parameters.sweep variable.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (numel(parameters.sweep)~=1)&&(numel(parameters.sweep)~=2)
    error('parameters.sweep array should have one or two elements.');
end
if ~isfield(parameters,'axis_units')
    error('axis units must be specified in parameters.axis_units variable.');
end
if ~ischar(parameters.axis_units)
    error('parameters.axis_units must be a character string.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)
    error('parameters.spins cell array must have exactly one element.');
end
end

% The only sin on earth is to do things badly.
%
% Ayn Rand

