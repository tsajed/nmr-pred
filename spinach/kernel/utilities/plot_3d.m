% Volume isosurface plotting utility with non-linear adaptive isosurface
% spacing.The function is useful for NMR data where small cross-peaks 
% must be adequately contoured next to large diagonal peaks. Syntax:
%
%     plot_3d(spin_system,spectrum,parameters,ncont,delta,k,signs)
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
%      for most 3D spectra is [0.02 0.2];
%      
%    * ncont is the number of contours, a reasonable value is 20.
%
%    * k controls the curvature of the contour spacing function: k=1 
%      corresponds to linear spacing and k>1 bends the spacing curve to
%      increase the sampling density near the baseline. A reasonable
%      value is 2;
%
%    * signs can be set to 'positive', 'negative' or 'both' - this will
%      cause the corresponding contours to be plotted.
%
% The following subfields are required inthe parameters structure:
%
%    parameters.sweep              three widths, Hz
%
%    parameters.spins              cell array with three character
%                                  strings specifying the working 
%                                  spins.
%
%    parameters.offset             three transmitter offsets, Hz
%
%    parameters.axis_units         axis units ('ppm','Hz','Gauss')
%
% i.kuprov@soton.ac.uk
% luke.edwards@ucl.ac.uk

function plot_3d(spin_system,spectrum,parameters,ncont,delta,k,signs)

% Check consistency
grumble(spectrum,parameters,ncont,delta,k,signs);

% Inform the user
report(spin_system,'plotting...');

% Determine data extents
xmax=max(spectrum(:)); xmin=min(spectrum(:)); 

% Compute contour levels
if (xmax>0)&&(strcmp(signs,'positive')||strcmp(signs,'both'))
    positive_contours=delta(2)*xmax*linspace(0,1,ncont).^k+xmax*delta(1);
else
    positive_contours=[];
end
if (xmin<0)&&(strcmp(signs,'negative')||strcmp(signs,'both'))
    negative_contours=delta(4)*xmin*linspace(0,1,ncont).^k+xmin*delta(3);
else
    negative_contours=[];
end
contours=[negative_contours(end:-1:1) positive_contours];

% Build axes and apply offsets
axis_f1=linspace(-parameters.sweep(1)/2,parameters.sweep(1)/2,size(spectrum,1))+parameters.offset(1);
axis_f2=linspace(-parameters.sweep(2)/2,parameters.sweep(2)/2,size(spectrum,2))+parameters.offset(2);
axis_f3=linspace(-parameters.sweep(3)/2,parameters.sweep(3)/2,size(spectrum,3))+parameters.offset(3);

% Convert the units
switch parameters.axis_units
    case 'ppm'
        axis_f1=1000000*(2*pi)*axis_f1/(spin(parameters.spins{1})*spin_system.inter.magnet);
        axis_f2=1000000*(2*pi)*axis_f2/(spin(parameters.spins{2})*spin_system.inter.magnet);
        axis_f3=1000000*(2*pi)*axis_f3/(spin(parameters.spins{3})*spin_system.inter.magnet);
        axis_label='ch. shift / ppm';
    case 'Gauss'
        axis_f1=10000*(spin_system.inter.magnet-2*pi*axis_f1/spin('E'));
        axis_f2=10000*(spin_system.inter.magnet-2*pi*axis_f2/spin('E'));
        axis_f3=10000*(spin_system.inter.magnet-2*pi*axis_f3/spin('E'));
        axis_label='magn. ind. / Gauss';
    case 'Hz'
        axis_f1=1*axis_f1+0; axis_f2=1*axis_f2+0; axis_f3=1*axis_f3+0;
        axis_label='lin. freq. / Hz';
    otherwise
        error('unknown axis units.');
end

% Plot the 3D spectrum
subplot(2,2,1,'replace'); hold('on');
[X,Y,Z]=meshgrid(axis_f2,axis_f1,axis_f3); 
for n=contours
    p=patch(isosurface(X,Y,Z,spectrum,n));
    set(p,'FaceColor','red','EdgeColor','none');
end
axis([min(axis_f2) max(axis_f2) min(axis_f1) max(axis_f1) min(axis_f3) max(axis_f3)]);
axis('square'); box('on'); grid('on'); camlight(); lighting('gouraud');
set(gca,'XDir','reverse','YDir','reverse','ZDir','reverse');
set(gca,'CameraPosition',[min(axis_f2)-max(axis_f2)...
                          min(axis_f1)-max(axis_f1)...
                          min(axis_f3)-max(axis_f3)]*100);
xlabel([parameters.spins{2} ' ' axis_label]); 
ylabel([parameters.spins{1} ' ' axis_label]); 
zlabel([parameters.spins{3} ' ' axis_label]); hold('off');

% Plot F3-F2 projection
subplot(2,2,3); proj_params=parameters;
proj_params.spins=proj_params.spins([3 2]);
proj_params.offset=proj_params.offset([3 2]);
proj_params.sweep=proj_params.sweep([3 2]);
proj_params.npoints=proj_params.npoints([3 2]);
proj_params.zerofill=proj_params.zerofill([3 2]);
plot_2d(spin_system,squeeze(sum(spectrum,1)),proj_params,20,delta,2,256,6,signs);
axis('square'); grid('on'); set(gca,'XDir','reverse','YDir','reverse');

% Plot F2-F1 projection
subplot(2,2,4); proj_params=parameters;
proj_params.spins=proj_params.spins([2 1]);
proj_params.offset=proj_params.offset([2 1]);
proj_params.sweep=proj_params.sweep([2 1]);
proj_params.npoints=proj_params.npoints([2 1]);
proj_params.zerofill=proj_params.zerofill([2 1]);
plot_2d(spin_system,squeeze(sum(spectrum,3)),proj_params,20,delta,2,256,6,signs);
axis('square'); grid('on'); set(gca,'XDir','reverse','YDir','reverse');

% Plot F3-F1 projection
subplot(2,2,2); proj_params=parameters;
proj_params.spins=proj_params.spins([3 1]);
proj_params.offset=proj_params.offset([3 1]);
proj_params.sweep=proj_params.sweep([3 1]);
proj_params.npoints=proj_params.npoints([3 1]);
proj_params.zerofill=proj_params.zerofill([3 1]);
plot_2d(spin_system,squeeze(sum(spectrum,2)),proj_params,20,delta,2,256,6,signs);
axis('square'); grid('on'); set(gca,'XDir','reverse','YDir','reverse');

% Make the figure bigger than the defalut
set(gcf,'Position',[100,100,1024,768]);

end

% Consistency enforcement
function grumble(spectrum,parameters,ncont,delta,k,signs)
if (~isnumeric(spectrum))||(~isreal(spectrum))||(ndims(spectrum)~=3)
    error('spectrum must be a real cube of numbers.');
end
if (~isfield(parameters,'offset'))
    error('offsets should be specified in parameters.offset variable.');
end
if numel(parameters.offset)~=3
    error('parameters.offset array should have three elements.');
end
if ~isfield(parameters,'sweep')
    error('sweep widths should be specified in parameters.sweep variable.');
end
if numel(parameters.sweep)~=3
    error('parameters.sweep array should have three elements.');
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
if ~iscell(parameters.spins)
    error('parameters.spins should be a cell array of character strings.');
end
if numel(parameters.spins)~=3
    error('parameters.spins cell array should have three elements.');
end
if (~isnumeric(ncont))||(~isscalar(ncont))||(~isreal(ncont))||(ncont<1)||(mod(ncont,1)~=0)
    error('ncont parameter must be a positive integer.');
end
if (~isnumeric(delta))||(numel(delta)~=4)||(~isreal(delta))||any(delta>1)||any(delta<0)
    error('delta parameter must be a vector with four real elements between 0 and 1.');
end
if (~isnumeric(k))||(~isscalar(k))||(~isreal(k))||(k<1)||(mod(k,1)~=0)
    error('k parameter must be a positive integer.');
end
if ~ischar(signs)
    error('signs parameter must be a character string.');
end
end

% The English are well known throughout the world for their lack of poli-
% tical scruples. They are experts at the art of hiding their misdeeds
% behind a facade of virtue. They have been at it for centuries, and it 
% has become such a part of their nature that they hardly notice it any 
% longer. They carry on with such a pious expression and deadly serious-
% ness that they even convince themselves that they are the exemplars of
% political virtue. They do not admit their hypocrisy to themselves. It 
% never happens that one Englishman says to another with a wink or a smi-
% le "We don't want to fool ourselves, do we now." They do not only beha-
% ve as if they were the model of piety and virtue - they really believe 
% that they are. That is both amusing and dangerous.
%
% Mickey Mouse

