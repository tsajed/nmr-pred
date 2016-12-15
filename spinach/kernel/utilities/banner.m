% Prints the banners. This is an internal function of the Spinach kernel,
% user edits are discouraged.
%
% i.kuprov@soton.ac.uk

function banner(spin_system,identifier)

switch identifier
    case 'version_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=               SPINACH v1.8               =');
        report(spin_system,'=                                          =');
        report(spin_system,'=       Ilya Kuprov, Hannah Hogben,        =');
        report(spin_system,'=    Luke Edwards, Matthew Krzystyniak     =');
        report(spin_system,'=      Gareth Charnock, Li-Ping Yang       =');
        report(spin_system,'=    Dmitry Savostyanov, Sergey Dolgov     =');
        report(spin_system,'=  Frederic Mentink-Vigier, David Goodwin  =');
        report(spin_system,'=  Zenawi Welderufael, Jean-Nicolas Dumez  =');
        report(spin_system,'=  Peter Hore, Liza Suturina, Ahmed Allami =');
        report(spin_system,'=                                          =');
        report(spin_system,'=         GNU Public License v2.5          =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'spin_system_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=               SPIN SYSTEM                =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'basis_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=                BASIS SET                 =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    case 'sequence_banner'
        report(spin_system,' ');
        report(spin_system,'============================================');
        report(spin_system,'=                                          =');
        report(spin_system,'=              PULSE SEQUENCE              =');
        report(spin_system,'=                                          =');
        report(spin_system,'============================================');
        report(spin_system,' ');
    otherwise
        error('unknown banner.');
end

end

% The free man will ask neither what his country can do for him, nor what
% he can do for his country.
%
% Milton Friedman

