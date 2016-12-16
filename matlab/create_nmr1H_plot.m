function [ spectrum ] = create_nmr1H_plot( input_args )
%% Create a function to plot 1H NMR given chemical shifts in ppm and J-coupling constants
%   Detailed explanation goes here
    % Spin system
    
sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H','1H','1H','1H','1H','1H',...
              '1H','1H','1H','1H'};

% Interactions
sys.magnet=5.9;
inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91...
                     3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72...
                     3.72 3.76 3.76 3.76};
inter.coupling.scalar{3,4}=12.1; 
inter.coupling.scalar{4,5}=3.1; 
inter.coupling.scalar{3,5}=1.0; 
inter.coupling.scalar{3,8}=1.0; 
inter.coupling.scalar{1,8}=1.0;
inter.coupling.scalar{6,7}=8.6; 
inter.coupling.scalar{5,8}=4.1; 
inter.coupling.scalar{7,9}=0.7; 
inter.coupling.scalar{7,10}=0.7; 
inter.coupling.scalar{9,10}=15.8;
inter.coupling.scalar{10,11}=9.8; 
inter.coupling.scalar{9,11}=8.1; 
inter.coupling.scalar{13,14}=1.5; 
inter.coupling.scalar{12,14}=0.9; 
inter.coupling.scalar{22,22}=0;

% Algorithmic options
sys.enable={'greedy'};
sys.disable={'colorbar'};
sys.tols.prox_cutoff=4.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-2';
bas.connectivity='scalar_couplings';
bas.space_level=1;
bas.sym_group={'S3','S3','S3'};
bas.sym_spins={[14 15 16],[17 18 19],[20 21 22]};

% Sequence parameters
parameters.angle=pi/2;
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

end

