% Contour plot slicing utility with non-linear adaptive contour spacing.
% Calls plot_3d() and allows slice extraction. Syntax:
%
%  plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)
%
% The following functions are used to compute contour levels:
%
%  cont_levs_pos=delta(2)*xmax*linspace(0,1,ncont).^k+xmax*delta(1);
%  cont_levs_neg=delta(2)*xmin*linspace(0,1,ncont).^k+xmin*delta(1);
%
% where:
%
%    * xmax and xmin are calculated from the spectrum;
%
%    * delta is the minimum and maximum elevation (as a fraction of total
%      intensity) of the contours above the baseline. A reasonable value 
%      for most 2D spectra is [0.02 0.2 0.02 0.2]. The first pair of num-
%      bers refers to the positive contours and the second pair to the
%      negative ones.
%
%    * ncont is the number of contours, a reasonable value is 20.
%
%    * k controls the curvature of the contour spacing function: k=1 
%      corresponds to linear spacing and k>1 bends the spacing curve to
%      increase the sampling density near the baseline. A reasonable
%      value is 2;
%
%    * ncol is a number of colors in the colormap (around 256 is fine);
%
%    * m is the curvature of the colormap: m=1 corresponds to a linear
%      color ramp into the red for positive contours and into the blue
%      for negative contours. A reasonable value for high-contrast
%      plotting is 6.
%
%    * signs can be set to 'positive', 'negative' or 'both' - this will
%      cause the corresponding contours to be plotted.
%
% The following subfields are required inthe parameters structure:
%
%    parameters.sweep              one or two sweep widths, Hz
%
%    parameters.spins              cell array with one ot two character
%                                  strings specifying the working spins.
%
%    parameters.offset             one or two transmitter offsets, Hz
%
%    parameters.axis_units         axis units ('ppm','Hz','Gauss')
%
% i.kuprov@soton.ac.uk

function slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)

% Do contour plotting
subplot(1,3,1); 
[f2,f1,S]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs);
title('2D spectrum'); drawnow();

% Switch off performance warning
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

% Get the spectrum extents
spec_min=min(spectrum(:));
spec_max=max(spectrum(:));

% Enter the slicing loop
while true()

    % Record pointer position
    subplot(1,3,1); [x,y]=ginput(1);
    
    % Compute axis grids
    [F1,F2]=ndgrid(f1,f2);
    
    % Create the interpolant
    S_int=griddedInterpolant(F1,F2,transpose(S),'spline');
    
    % Compute the traces
    trace_f1=S_int(ones(size(f2))*x,f2);
    trace_f2=S_int(f1,ones(size(f1))*y);
    
    % Call the 1D plotting routine for F1
    parameters_f1=parameters;
    parameters_f1.spins=parameters.spins(1);
    parameters_f1.offset=parameters.offset(1);
    parameters_f1.sweep=parameters.sweep(1);
    parameters_f1.zerofill=parameters.zerofill(1);
    subplot(1,3,2); plot_1d(spin_system,trace_f1,parameters_f1);
    set(gca,'YLim',[spec_min spec_max]); title('F1 slice'); drawnow();
    
    
    % Call the 1D plotting routine for F2
    parameters_f2=parameters;
    if numel(parameters.spins)==2
        parameters_f2.spins=parameters.spins(2);
    end
    if numel(parameters.offset)==2
        parameters_f2.offset=parameters.offset(2);
    end
    parameters_f2.sweep=parameters.sweep(2);
    parameters_f2.zerofill=parameters.zerofill(2);
    subplot(1,3,3); plot_1d(spin_system,trace_f2,parameters_f2);
    set(gca,'YLim',[spec_min spec_max]); title('F2 slice'); drawnow();
    
end

end

% There's a great deal of difference between an eager man who wants
% to read a book and a tired man who wants a book to read.
%
% Gilbert K. Chesterton

