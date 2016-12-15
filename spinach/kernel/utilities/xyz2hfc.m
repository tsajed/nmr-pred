% Converts point electron and nuclear coordinates (Angstroms) 
% into a hyperfine interaction tensor in Hz. Syntax:
%
%                A=xyz2hfc(e_xyz,n_xyz,isotope)
%
% i.kuprov@soton.ac.uk

function A=xyz2hfc(e_xyz,n_xyz,isotope)

% Fundamental constants
hbar=1.054571800e-34; mu0=4*pi*1e-7;

% Get magnetogyric ratios
gamma_e=spin('E'); gamma_n=spin(isotope);

% Get the distance
distance=norm(e_xyz-n_xyz);
        
% Get the ort
ort=(e_xyz-n_xyz)/distance;
        
% Compute the dipolar interaction prefactor
D=gamma_e*gamma_n*hbar*mu0/(4*pi*(distance*1e-10)^3);
        
% Compute the dipolar coupling matrix
A=D*[1-3*ort(1)*ort(1)   -3*ort(1)*ort(2)   -3*ort(1)*ort(3);
      -3*ort(2)*ort(1)  1-3*ort(2)*ort(2)   -3*ort(2)*ort(3);
      -3*ort(3)*ort(1)   -3*ort(3)*ort(2)  1-3*ort(3)*ort(3)];
        
end

% "English people as a whole have a rooted distrust of total
%  abstainers as politicians."
%
% The Very Rev'd H. Hensley Henson, then Dean of Durham, 
% in a letter to the Times, 4th January 1916.

