% Fits a Karplus curve to a Gaussian dihedral angle scan.
%
% <http://spindynamics.org/wiki/index.php?title=Karplus_fit.m>

function [A,B,C]=karplus_fit(dir_path,atoms)

% Get all log files in the directory
logfiles=dir([dir_path '/*.log']);

% Get the arrays going
phi=[]; J=[]; E=[];

% Extract parameters
for n=1:numel(logfiles)
    try
        props=gparse([dir_path '/' logfiles(n).name]);
        current_phi=dihedral(props.std_geom(atoms(1),:),props.std_geom(atoms(2),:),...
                             props.std_geom(atoms(3),:),props.std_geom(atoms(4),:));
        current_j=props.j_couplings(atoms(1),atoms(4));
        phi=[phi current_phi]; J=[J current_j]; E=[E props.energy]; %#ok<AGROW>
    catch
        disp(['Gaussian log file ' logfiles(n).name ' could not be interpreted, skipped.']);
    end
end

% Rotate phi into [0,360]
phi=mod(phi,360);

% Compute basis cosines
cosine_2=cosd(phi).^2;
cosine_1=cosd(phi);
cosine_0=ones(size(phi));

% Run linear least squares
result=[cosine_2' cosine_1' cosine_0']\J';

% Write the answer
A=result(1); B=result(2); C=result(3);

% Plot Karplus curve
figure(); subplot(2,1,1); plot(phi,J,'ro');
psi=linspace(0,360,128); hold on; 
plot(psi,A*cosd(psi).^2+B*cosd(psi)+C,'b-');
xlabel('Dihedral angle, degrees');
ylabel('J-coupling, Hz'); axis tight;

% Plot probability histogram
E=hartree2joule(E); E=E-min(E(:));
P=exp(-E/(8.31*298)); P=P/sum(P);
[phi,order]=sort(phi); subplot(2,1,2);
stairs(phi,P(order));
xlabel('Dihedral angle, degrees');
ylabel('Probability at 298 K'); axis tight;

end

% If men learn this, it will implant forgetfulness in their souls; they will
% cease to exercise memory because they rely on that which is written [...]
% it is no true wisdom that you offer your disciples, but only its semblance,
% for by telling them of many things without teaching them you will make them
% seem to know much, while for the most part they know nothing, and as men 
% filled, not with wisdom but with the conceit of wisdom, they will be a bur-
% den to their fellows.
%
% Plato (circa 429-347 BCE), about reading and writing.

