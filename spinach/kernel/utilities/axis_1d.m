% Generates axes for plotting.
%
% i.kuprov@soton.ac.uk

function [ax,ax_label]=axis_1d(spin_system,parameters)

% Build the axis and apply the offset
if numel(parameters.sweep)==1
    ax=linspace(-parameters.sweep/2,parameters.sweep/2,parameters.zerofill)+parameters.offset;
else
    ax=linspace(parameters.sweep(1),parameters.sweep(2),parameters.zerofill);
end

% Convert the units if necessary
switch parameters.axis_units
    case 'ppm'
        ax_label=[parameters.spins{1} ' chemical shift / ppm'];
        ax=1000000*(2*pi)*ax/(spin(parameters.spins{1})*spin_system.inter.magnet);
    case 'Gauss'
        ax_label='Magnetic induction / Gauss';
        ax=10000*(spin_system.inter.magnet-2*pi*ax/spin('E'));
    case 'Gauss-labframe'
        ax_label='Magnetic induction / Gauss';
        ax=10000*(2*pi*ax/spin('E'));
    case 'mT'
        ax_label='Magnetic induction / mT';
        ax=1000*(spin_system.inter.magnet-2*pi*ax/spin('E'));
    case 'mT-labframe'
        ax_label='Magnetic induction / mT';
        ax=1000*(2*pi*ax/spin('E'));
    case 'T'
        ax_label='Magnetic induction / Tesla';
        ax=spin_system.inter.magnet-2*pi*ax/spin('E');
    case 'Hz'
        ax_label=[parameters.spins{1} ' offset linear frequency / Hz'];
        ax=1*ax;    
    case 'kHz'
        ax_label=[parameters.spins{1} ' offset linear frequency / kHz'];
        ax=0.001*ax;
    case 'MHz'
        ax_label=[parameters.spins{1} ' offset linear frequency / MHz'];
        ax=0.000001*ax;
    case 'MHz-labframe'
        ax_label=[parameters.spins{1} ' linear frequency / MHz'];
        ax=0.000001*(ax-spin_system.inter.magnet*spin(parameters.spins{1})/(2*pi));
    case 'GHz'
        ax_label=[parameters.spins{1} ' offset linear frequency / GHz'];
        ax=0.000000001*ax;
    case 'GHz-labframe'
        ax_label=[parameters.spins{1} ' linear frequency / GHz'];
        ax=0.000000001*(ax-spin_system.inter.magnet*spin(parameters.spins{1})/(2*pi));
    otherwise
        error('unknown axis units.');
end

end

% The God of the Old Testament is arguably the most unpleasant character in 
% all fiction: jealous and proud of it; a petty, unjust, unforgiving control
% freak; a vindictive, bloodthirsty ethnic cleanser; a misogynistic, homopho-
% bic, racist, infanticidal, genocidal, filicidal, pestilential, megalomania-
% cal, sadomasochistic, capriciously malevolent bully.
%
% Richard Dawkins, "The God Delusion"

