% 2D spectral integrator. Calls plot_2d.m and then launches either an
% interactive integration procedure (if no range file name is given) 
% and records the range data into a file, or runs an automatic integ-
% ration if the range file name from a previous run is supplied.
% Syntax:
% 
%    int_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,...
%           m,signs,filename)
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

function int_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs,filename)

% Do contour plotting
[f2,f1,S]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs);

% Switch off performance warning
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

% Check if the file exists
if ~exist('filename','var')
    
    % Proceed with interactive integration
    while true()
        
        % Report to the user
        disp('Interactive integration, define the box by clicking on its opposite corners.');
        
        % Record pointer position
        [corners_x,corners_y]=ginput(2);
        
        % Compute axis grids
        [F1,F2]=ndgrid(f1,f2);
        
        % Create the interpolant
        S_int=griddedInterpolant(F1,F2,transpose(S),'spline');
        
        % Make function handle
        fhandle=@(x,y)S_int(x,y);
        
        % Compute the integral
        I=integral2(fhandle,corners_x(1),corners_x(2),corners_y(1),corners_y(2),'RelTol',1e-3,'AbsTol',1e-3);
        
        % Report to the user
        disp(['X range: ' num2str(min(corners_x)) ' to ' num2str(max(corners_x))]);
        disp(['Y range: ' num2str(min(corners_y)) ' to ' num2str(max(corners_y))]);
        disp(['Integral: ' num2str(abs(I))]);
        
    end
    
elseif exist(filename,'file')
    
    % Load ranges from file
    disp('Found the ranges file, integrals:'); load(filename,'ranges');
    
    % Perform automatic integration
    for n=1:size(ranges,1) %#ok<NODEF>
        
        % Get pointer position
        corners_x=ranges{n,1}; corners_y=ranges{n,2};
        
        % Compute axis grids
        [F1,F2]=ndgrid(f1,f2);
        
        % Create the interpolant
        S_int=griddedInterpolant(F1,F2,transpose(S),'spline');
        
        % Make function handle
        fhandle=@(x,y)S_int(x,y);
        
        % Compute the integral
        I=integral2(fhandle,corners_x(1),corners_x(2),corners_y(1),corners_y(2),'RelTol',1e-3,'AbsTol',1e-3);
        
        % Report to the user
        disp(I);
        
    end
    
else
    
    % Record ranges into a file
    n=1;
    while true()
        
        % Get mouse input
        [ranges{n,1},ranges{n,2}]=ginput(2); %#ok<AGROW>
        
        % Write the file
        save(filename,'ranges');
        
        % Give feedback
        disp(['Recorded X range: ' num2str(min(ranges{n,1})) ' to ' num2str(max(ranges{n,1}))]);
        disp(['Recorded Y range: ' num2str(min(ranges{n,2})) ' to ' num2str(max(ranges{n,2}))]);
        
        % Increment counter
        n=n+1;
        
    end
    
end

end

% Beware of the person of one book.
%
% Thomas Aquinas

