% Shaped pulse function using Cartesian coordinates.
%
% <http://spindynamics.org/wiki/index.php?title=Shaped_pulse_xy.m>

function [rho,P]=shaped_pulse_xy(spin_system,drift,controls,amplitudes,time_grid,rho)

% Check consistency
grumble(drift,controls,amplitudes,time_grid,rho);

% Decide the methods
switch nargout

    case 1
        
        % Run Krylov propagation
        for n=1:numel(time_grid)
            
            % Generate the evolution step operator
            slice_operator=drift;
            for k=1:numel(controls)
                slice_operator=slice_operator+amplitudes{k}(n)*controls{k};
            end
            
            % Apply the evolution slice
            rho=step(spin_system,slice_operator,rho,time_grid(n));
            
        end
        
    case 2
        
        % Get the propagator going
        P=speye(size(drift));
        
        % Compute the pulse propagator
        parfor n=1:numel(time_grid)
            
            % Generate the evolution step operator
            slice_operator=drift;
            for k=1:numel(controls)
                slice_operator=slice_operator+amplitudes{k}(n)*controls{k}; %#ok<PFBNS>
            end
            
            % Apply the evolution slice
            P=propagator(spin_system,slice_operator,time_grid(n))*P;
            
        end
        
        % Apply pulse propagator
        rho=P*rho;

end

end

% Consistency enforcement
function grumble(drift,controls,amplitudes,time_grid,rho) %#ok<INUSD>

% to be written

end

% When IK proposed the fibre etching technique featured in the 2004
% JMR paper (http://dx.doi.org/10.1016/j.jmr.2004.08.017), he could
% see terror in his supervisors's eyes - Peter Hore reasonably tho-
% ught that the notoriously eccentric Russian student could not pos-
% sibly be trusted with boiling hydrofluoric acid in a high-power 
% laser lab. The Chemistry Department safety officer held a similar
% view. It is not entirely clear how IK got hold of several millili-
% ters of concentrated HF and a heat gun on a Saturday night in the
% PTCL Teaching Lab, but the photographs of the resulting fibre tip
% were left on Peter's table on Monday. The paper was accepted by 
% JMR without revisions.

