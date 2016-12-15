% Returns a rotor stack of Liouvillians or Hamiltonians. The stack is 
% needed for the traditional style calculation of MAS dynamics. Syntax:
%
%            L=rotor_stack(spin_system,parameters,assumptions)
%
% Parameters:
%
%   parameters.axis     - spinning axis, given as a normalized
%                         3-element vector
%
%   parameters.offset   - a cell array giving transmitter off-
%                         sets in Hz on each of the spins listed
%                         in parameters.spins array
%
%   parameters.max_rank - maximum harmonic rank to retain in
%                         the solution (increase till conver-
%                         gence is achieved, approximately
%                         equal to the number of spinning si-
%                         debands in the spectrum)
%
%   parameters.rframes  - rotating frame specification, e.g.
%                         {{'13C',2},{'14N,3}} requests second
%                         order rotating frame transformation
%                         with respect to carbon-13 and third
%                         order rotating frame transformation
%                         with respect to nitrogen-14. When
%                         this option is used, the assumptions
%                         on the respective spins should be
%                         laboratory frame.
%
%   parameters.orientation - the orientation of the spin system at rotor
%                            phase zero, a vector of three Euler angles
%                            in radians.
%
% Note: relaxation and chemical kinetics are not included.
%
% i.kuprov@soton.ac.uk

function L=rotor_stack(spin_system,parameters,assumptions)

% Get the Hamiltonian
[H,Q]=hamiltonian(assume(spin_system,assumptions));

% Apply offsets
H=frqoffset(spin_system,H,parameters);

% Get rotor axis orientation
[phi,theta,~]=cart2sph(parameters.axis(1),parameters.axis(2),parameters.axis(3)); 
theta=pi/2-theta; D_lab2rot=euler2wigner(phi,theta,0); D_rot2lab=D_lab2rot';

% Compute rotor angles
[rotor_angles,~]=fourdif(2*parameters.max_rank+1,1);

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Set crystallite orientation
D_crystal=euler2wigner(parameters.orientation);

% Preallocate Liouvillian blocks
L=cell(2*parameters.max_rank+1,1);

% Build Liouvillian blocks
parfor n=1:(2*parameters.max_rank+1) %#ok<*PFBNS>
    
    % Get the rotation
    D_rotor=euler2wigner(0,0,rotor_angles(n));
    if strcmp(parameters.masframe,'magnet')
        D=D_rot2lab*D_rotor*D_lab2rot*D_crystal;
    elseif strcmp(parameters.masframe,'rotor')
        D=D_rot2lab*D_rotor*D_crystal;
    else
        D=0; error('unknown MAS frame.'); %#ok<NASGU>
    end
        
    % Build the block
    L{n}=H;
    for k=1:5
        for m=1:5
            L{n}=L{n}+D(k,m)*Q{k,m};
        end
    end
    L{n}=(L{n}+L{n}')/2;
    
    % Apply the rotating frame
    for k=1:numel(parameters.rframes) 
        L{n}=rotframe(spin_system,C{k},L{n},parameters.rframes{k}{1},parameters.rframes{k}{2});
    end
    
    % Clean up the result
    L{n}=clean_up(spin_system,L{n},spin_system.tols.liouv_zero);
    
end

end

% Never complain and never explain. 
%
% Benjamin Disraeli

