% A basic hydrodynamics infrastructure provider.
%
% i.kuprov@soton.ac.uk
% a.j.allami@soton.ac.uk

function [Fx,Fy,Fz]=hydrodynamics(parameters)

% Build derivative operators
switch parameters.cond{1}
    
    case 'period'
        
        % Finite-difference derivatives
        if numel(parameters.npts)==1
            Dx=fdmat(parameters.npts(1),parameters.cond{2},1)/(parameters.dims(1)/parameters.npts(1));
        end
        if numel(parameters.npts)==2
            Dx=fdmat(parameters.npts(1),parameters.cond{2},1)/(parameters.dims(1)/parameters.npts(1));
            Dy=fdmat(parameters.npts(2),parameters.cond{2},1)/(parameters.dims(2)/parameters.npts(2));
        end
        if numel(parameters.npts)==3
            Dx=fdmat(parameters.npts(1),parameters.cond{2},1)/(parameters.dims(1)/parameters.npts(1));
            Dy=fdmat(parameters.npts(2),parameters.cond{2},1)/(parameters.dims(2)/parameters.npts(2));
            Dz=fdmat(parameters.npts(3),parameters.cond{2},1)/(parameters.dims(3)/parameters.npts(3));
        end
        
   case 'fourier'
        
        % Fourier derivatives
        if numel(parameters.npts)==1
            [~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx;
        end
        if numel(parameters.npts)==2
            [~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx;
            [~,Dy]=fourdif(parameters.npts(2),1); Dy=(2*pi/parameters.dims(2))*Dy;
        end
        if numel(parameters.npts)==3
            [~,Dx]=fourdif(parameters.npts(1),1); Dx=(2*pi/parameters.dims(1))*Dx;
            [~,Dy]=fourdif(parameters.npts(2),1); Dy=(2*pi/parameters.dims(2))*Dy;
            [~,Dz]=fourdif(parameters.npts(3),1); Dz=(2*pi/parameters.dims(3))*Dz;
        end
        
    otherwise
        
        % Complain and bomb out
        error('unrecognized derivative operator type.');
        
end

% Kron up derivative operators
if numel(parameters.npts)==1
    Fx=-1i*Dx; Fy=[]; Fz=[];
end
if numel(parameters.npts)==2
    Fx=-1i*kron(kron(speye(parameters.npts(2)),Dx),speye(size(1)));
    Fy=-1i*kron(kron(Dy,speye(parameters.npts(1))),speye(size(1))); Fz=[];
end
if numel(parameters.npts)==3
    Fx=-1i*kron(kron(kron(speye(parameters.npts(3)),speye(parameters.npts(2))),Dx),speye(size(1)));
    Fy=-1i*kron(kron(kron(speye(parameters.npts(3)),Dy),speye(parameters.npts(1))),speye(size(1)));
    Fz=-1i*kron(kron(kron(Dz,speye(parameters.npts(2))),speye(parameters.npts(1))),speye(size(1)));
end

end

% As for butter versus margarine, I trust cows
% more than chemists.
%
% Joan Gussow

