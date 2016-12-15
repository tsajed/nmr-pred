% Rotating frame transformation with respect to a specified
% group of spins to specified order in perturbation theory.
% Syntax:
%
%          H=rotframe(spin_system,H0,H,isotope,order)
%
% Parameters:
%
%    H0     - carrier Hamiltonian with respect to which the
%             rotating frame transformation is to be done
%
%    H      - laboratory frame Hamiltonian H0+H1 that is to
%             be transformed into the rotating frame 
%
%    isotope - string, such as '1H', specifying the spins
%              with respect to which the transformation is
%              being computed
%
%    order   - perturbation theory order in the rotating
%              frame transformation
%
% i.kuprov@soton.ac.uk

function Hr=rotframe(spin_system,H0,H,isotope,order)

% Check consistency
grumble(H0,H,order);

% Compute the period
switch spin_system.bas.formalism
    
    case {'zeeman-liouv','sphten-liouv'}
        
        % Liouville space period for H0
        T=-2*pi/(spin(isotope)*spin_system.inter.magnet);
        
    case {'zeeman-hilb'}
        
        % Hilbert space period for H0
        T=-4*pi/(spin(isotope)*spin_system.inter.magnet);
        
end

% Run the interaction representation transformation
Hr=intrep(spin_system,H0,H,T,order);

end

% Consistency enforcement
function grumble(H0,H,order)
if (~ishermitian(H))||(~ishermitian(H0))
    error('both H and C must be Hermitian.');
end
if ((~isreal(order))||(order<1)||(mod(order,1)~=0))&&(~isinf(order))
    error('unsupported rotating frame transformation theory order.');
end
end

% Finally, I acknowledge the financial support of EPSRC in its best,
% that is, the responsive mode. This Nobel Prize would be absolutely
% impossible without this mode.  [...]  However, I can offer no nice
% words for the EU Framework Programmes which, except for the Europe-
% an Research Council, can be praised only by europhobes for discre-
% diting the whole idea of an effectively working Europe.
%
% Andre Geim's Nobel Lecture, 2010

