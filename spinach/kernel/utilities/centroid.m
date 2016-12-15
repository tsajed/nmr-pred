% Finds the centroid of a 3D probability density in a cube.
%
% i.kuprov@outlook.com

function [x,y,z]=centroid(probden,ranges)

% Get coordinate arrays
[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(probden,1)),...
               linspace(ranges(3),ranges(4),size(probden,2)),...
               linspace(ranges(5),ranges(6),size(probden,3)));

% Get the norm
n=sum(sum(sum(probden)));

% Get centroid coordinates
x=sum(sum(sum(X.*probden)))/n;
y=sum(sum(sum(Y.*probden)))/n;
z=sum(sum(sum(Z.*probden)))/n;

end

% Life is warfare.
%
% Lucius Annaeus Seneca

