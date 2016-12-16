%function [ spectrum ] = create_nmr1H_plot( sys, inter )
%% Create plot 1H NMR given chemical shifts in ppm and J-coupling constants
%   Detailed explanation goes here
    % Spin system
    
% sys.isotopes={'1H','1H','1H','1H','1H'};

% Interactions
sys.magnet=5.9;
% inter.zeeman.scalar={2.0, 2.0, 3.0, 3.0, 3.0};

% inter.coupling.scalar{1,2}=1; 
% inter.coupling.scalar{1,3}=1; 
% inter.coupling.scalar{1,4}=1; 
% inter.coupling.scalar{1,5}=1; 
% inter.coupling.scalar{2,3}=1; 
% inter.coupling.scalar{2,4}=1; 
% inter.coupling.scalar{2,5}=1; 
% inter.coupling.scalar{3,4}=1; 
% inter.coupling.scalar{4,5}=1; 
% inter.coupling.scalar{5,5}=0; 

% Algorithmic options
sys.enable={'greedy'};
sys.disable={'colorbar'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;
% bas.sym_group={'S3','S3','S3'};
% bas.sym_spins={[14 15 16],[17 18 19],[20 21 22]};

% Sequence parameters
% parameters.angle=pi/2;
parameters.offset=1200;
parameters.sweep=2000;
parameters.npoints=512;
parameters.zerofill=2048;
parameters.spins={'1H'};
parameters.axis_units='ppm';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

parameters.rho0=state(spin_system,'L+','1H','cheap');
parameters.coil=state(spin_system,'L+','1H','cheap');

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodization
fid=apodization(fid,'crisp-1d');

% Fourier transform
spectrum=fftshift(fft(fid,parameters.zerofill));

% Destreaking
%spectrum=destreak(spectrum);

% Plotting
figure(); 
%plot_2d(spin_system,real(spectrum),parameters,20,[0.01 0.1 0.01 0.1],2,256,6,'both');

plot_1d(spin_system,real(spectrum),parameters);

%end

