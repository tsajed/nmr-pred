% ORCA cube file parser. Extracts the normalised probability density and
% the associated metric information from ORCA cube files. Syntax:
%
%                 [density,ext,dx,dy,dz]=ocparse(filename)
%
% Outputs:
%
%    density  - probability density cube with dimensions
%               ordered as [X Y Z]
%
%    ext      - grid extents in Angstrom, ordered as
%               [xmin xmax ymin ymax zmin zmax]
%
%    dx,dy,dz - grid steps in the three directions, Angstrom
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk

function [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)

% Inform the user
disp('Parsing ORCA cube...');

% Parse the file
A=importdata(filename);

% Get the number of points along [x y z]
npts=str2num(A.textdata{2}); %#ok<ST2NM>

% Make the probability density cube
density=reshape(abs(A.data),npts(2),npts(3),npts(1));
density=permute(density,[2 3 1]);

% Get the corner coordinates 
corner_xyz=str2num(A.textdata{3}); %#ok<ST2NM>

% Get the grid spacing vector
dxdydz=str2num(A.textdata{4}); %#ok<ST2NM>
dx=dxdydz(1); dy=dxdydz(2); dz=dxdydz(3);

% Integrate and normalise probability density
total_prob=simps(simps(simps(density)))*dx*dy*dz;
density=density/total_prob;

% Compute grid extents
ext=[corner_xyz(1) (corner_xyz(1)+(npts(1)-1)*dx)...
     corner_xyz(2) (corner_xyz(2)+(npts(2)-1)*dy)...
     corner_xyz(3) (corner_xyz(3)+(npts(3)-1)*dz)];
 
% Pad the density with zeros
density=padarray(density,pad_factor*size(density),0,'both');
ext=(2*pad_factor+1)*ext;

end

% There is a cult of ignorance in the United States, and there
% has always been. The strain of anti-intellectualism has been
% a constant thread winding its way through our political and
% cultural life, nurtured by the false notion that democracy
% means that 'my ignorance is just as good as your knowledge'. 
%
% Isaac Asimov

