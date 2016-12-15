% Finds the maximum of the probability density.
%
% i.kuprov@outlook.com

function [x,y,z]=probmax(probden,ranges)

% Get coordinate arrays
[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(probden,1)),...
               linspace(ranges(3),ranges(4),size(probden,2)),...
               linspace(ranges(5),ranges(6),size(probden,3)));

% Get the max
[~,index]=max(probden(:));

% Get maximum coordinates
x=X(index); y=Y(index); z=Z(index);

end

% What exactly is your "fair share" of what 
% someone else has worked for?
%
% Thomas Sowell

