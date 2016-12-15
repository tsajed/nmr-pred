% Spherical grid direct product. Tiles one grid using the rotations of
% the other. Grids should be supplied using Euler angles in three col-
% umns [alphas betas gammas] in radians. Syntax:
%
%    [angles,weights]=grid_kron(angles1,weights1,angles2,weights2)
%
% Parameters:
%
%    angles1 - angles of the first grid, as [alpha beta gamma], rad
%
%   weights1 - weights of the first grid
%
%    angles2 - angles of the second grid, as [alpha beta gamma], rad
%
%   weights1 - weights of the second grid
%
% Outputs:
%
%     angles - angles of the product grid, as [alpha beta gamma], rad
%
%    weights - weights of the product grid
%
% i.kuprov@soton.ac.uk

function [angles,weights]=grid_kron(angles1,weights1,angles2,weights2)

% Convert both grids to quaternions
quats1=angle2quat(angles1(:,1),angles1(:,2),angles1(:,3),'ZYZ');
quats2=angle2quat(angles2(:,1),angles2(:,2),angles2(:,3),'ZYZ');

% Build a table of quaternion products
quats=[kron(quats1,ones(size(quats2,1),1)) kron(ones(size(quats1,1)),quats2)];

% Multiply up quaternions
quats=[quats(:,1).*quats(:,5)-quats(:,2).*quats(:,6)-quats(:,3).*quats(:,7)-quats(:,4).*quats(:,8),...
       quats(:,1).*quats(:,6)+quats(:,2).*quats(:,5)+quats(:,3).*quats(:,8)-quats(:,4).*quats(:,7),...
       quats(:,1).*quats(:,7)-quats(:,2).*quats(:,8)+quats(:,3).*quats(:,5)+quats(:,4).*quats(:,6),...
       quats(:,1).*quats(:,8)+quats(:,2).*quats(:,7)-quats(:,3).*quats(:,6)+quats(:,4).*quats(:,5)];
   
% Convert quaternions into angles
[alphas,betas,gammas]=quat2angle(quats,'ZYZ'); angles=real([alphas betas gammas]);

% Tile the weights
weights=kron(weights1,weights2);

end

% Dostoevsky's lack of taste, his monotonous dealings with persons
% suffering with pre-Freudian complexes, the way he has of wallow-
% ing in the tragic misadventures of human dignity -- all this is
% difficult to admire.
%
% Vladimir Nabokov

