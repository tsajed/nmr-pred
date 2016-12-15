% Crops 2D spectra to user-specified ranges (in ppm).
%
% i.kuprov@soton.ac.uk

function [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)

% Accommodate homonuclear 2D sequences
if numel(parameters.spins)==1,  parameters.spins=[parameters.spins parameters.spins]; end
if numel(parameters.offset)==1, parameters.offset=[parameters.offset parameters.offset]; end
if numel(parameters.sweep)==1,  parameters.sweep=[parameters.sweep parameters.sweep]; end

% Build axes and apply offsets
axis_f1_hz=linspace(-parameters.sweep(1)/2,parameters.sweep(1)/2,size(spec,1))+parameters.offset(1);
axis_f2_hz=linspace(-parameters.sweep(2)/2,parameters.sweep(2)/2,size(spec,2))+parameters.offset(2);

% Convert the units 
axis_f1_ppm=1000000*(2*pi)*axis_f1_hz/(spin(parameters.spins{1})*spin_system.inter.magnet);
axis_f2_ppm=1000000*(2*pi)*axis_f2_hz/(spin(parameters.spins{2})*spin_system.inter.magnet);

% Find array bounds
l_bound_f1=find(axis_f1_ppm>crop_ranges{1}(1),1);
r_bound_f1=find(axis_f1_ppm>crop_ranges{1}(2),1);
l_bound_f2=find(axis_f2_ppm>crop_ranges{2}(1),1);
r_bound_f2=find(axis_f2_ppm>crop_ranges{2}(2),1);

% Find the new offsets
parameters.offset=[(axis_f1_hz(l_bound_f1)+axis_f1_hz(r_bound_f1))/2 ...
                   (axis_f2_hz(l_bound_f2)+axis_f2_hz(r_bound_f2))/2];

% Find the new sweeps
parameters.sweep= [(axis_f1_hz(r_bound_f1)-axis_f1_hz(l_bound_f1)) ...
                   (axis_f2_hz(r_bound_f2)-axis_f2_hz(l_bound_f2))];
               
% Update the point counts
parameters.zerofill=[(r_bound_f1-l_bound_f1) (r_bound_f2-l_bound_f2)];

% Cut the spectrum
spec=spec(l_bound_f1:r_bound_f1,l_bound_f2:r_bound_f2);

end

% There is long term price instability for lanthanides, as most come 
% from China, but are produced with little regard for the environment
% or workers' health. Attempts to improve conditions led to a 3000%
% price increase for dysprosium in 2011.
%
% EPSRC reviewer, on one of IK's grant applications

