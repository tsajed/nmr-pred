% Translates a stationary 3D velocity field into a Fokker-Planck
% evolution generator.
%
% a.j.allami@soton.ac.uk
% i.kuprov@soton.ac.uk

function F=v2fplanck(U,V,W,parameters)

% Get the translation generators
[Fx,Fy,Fz]=hydrodynamics(parameters);

% Build the Fokker-Planck flow generator
F=spdiags(Fx*U(:)+Fy*V(:)+Fz*W(:),0,prod(parameters.npts),prod(parameters.npts))+...
  spdiags(U(:),0,prod(parameters.npts),prod(parameters.npts))*Fx+...
  spdiags(V(:),0,prod(parameters.npts),prod(parameters.npts))*Fy+...
  spdiags(W(:),0,prod(parameters.npts),prod(parameters.npts))*Fz;

% Add the diffusion term
F=F-1i*parameters.diffc*(Fx*Fx+Fy*Fy+Fz*Fz);

end

% "He's spending a year dead for tax reasons."
%
% Douglas Adams, The Hitchhiker's Guide to the Galaxy

