% Reads SIMPSON spin system specification file and converts the
% information into Spinach data structures. Syntax:
%
%                 [sys,inter]=s2spinach(filename)
%
% i.kuprov@soton.ac.uk

function [sys,inter]=s2spinach(filename)

% Check consistency
grumble(filename);

% Read the file
fid=fopen(filename);
simpson_file=textscan(fid,'%s','delimiter','\n');
simpson_file=simpson_file{:}; fclose(fid);

% Locate the boundaries of the spinsys{} field
for n=1:length(simpson_file)
    current_line=deblank(char(simpson_file{n}));
    if (length(current_line)>=7)&&strcmp(current_line(1:7),'spinsys')
        spinsys_start=n;
        for k=(n+1):length(simpson_file)
            current_line=deblank(char(simpson_file{k}));
            if strcmp(current_line,'}')
                spinsys_end=k; break
            end
        end
    end
end
 
% Tell the user if we found the spinsys{} section
if exist('spinsys_start','var')&&exist('spinsys_end','var')
    disp(['SIMPSON spinsys{} section located between lines ' num2str(spinsys_start) ' and ' num2str(spinsys_end) '.']);
else
    error('could not locate SIMPSON spinsys{} section');
end

% Search for the isotope information
for n=spinsys_start:spinsys_end
    current_line=deblank(char(simpson_file{n}));
    if (length(current_line)>=6)&&strcmp(current_line(1:6),'nuclei')
        sys.isotopes=textscan(current_line(7:end),'%s');
        sys.isotopes=sys.isotopes{:}';
    end
end

% Preallocate the coupling arrays
nspins=length(sys.isotopes);
inter.zeeman.matrix=mat2cell(zeros(3*nspins,3),3*ones(nspins,1))';
inter.coupling.matrix=mat2cell(zeros(3*nspins),3*ones(nspins,1),3*ones(nspins,1));

% Absorb the couplings
for n=spinsys_start:spinsys_end
    
    % Kill the spaces at either end
    current_line=deblank(char(simpson_file{n}));
    
    % Parse the specification
    if (length(current_line)>=5)&&strcmp(current_line(1:5),'shift')
        
        % Bomb out if the Zeeman tensor is not in ppm
        if nnz(current_line(6:end)=='p')~=2
            error('only ppm specification (# #p #p # # # #) is supported for Zeeman tensors.');
        end
        
        % Get the ppm Zeeman tensor information
        zeeman_line=textscan(current_line(6:end),'%n %np %np %n %n %n %n');
        eigvals=[zeeman_line{2}-zeeman_line{3}*(1+zeeman_line{4})/2
                 zeeman_line{2}-zeeman_line{3}*(1-zeeman_line{4})/2
                 zeeman_line{3}+zeeman_line{2}];
        S=euler2dcm(pi*[zeeman_line{5} zeeman_line{6} zeeman_line{7}]/180);
        inter.zeeman.matrix{zeeman_line{1}}=inter.zeeman.matrix{zeeman_line{1}}+S*diag(eigvals)*S';
        
    elseif (length(current_line)>=6)&&strcmp(current_line(1:6),'dipole')
        
        % Get the dipole-dipole coupling information
        dipole_line=textscan(current_line(7:end),'%n %n %n %n %n %n');
        eigvals=[-dipole_line{3}/2 -dipole_line{3}/2 dipole_line{3}];
        S=euler2dcm(pi*[dipole_line{4} dipole_line{5} dipole_line{6}]/180);
        inter.coupling.matrix{dipole_line{1},dipole_line{2}}=inter.coupling.matrix{dipole_line{1},dipole_line{2}}+S*diag(eigvals)*S';
        
    elseif (length(current_line)>=9)&&strcmp(current_line(1:9),'jcoupling')
        
        % Get the generic coupling information
        coupling_line=textscan(current_line(10:end),'%n %n %n %n %n %n %n %n');
        eigvals=[coupling_line{3}-coupling_line{4}*(1+coupling_line{5})/2
                 coupling_line{3}-coupling_line{4}*(1-coupling_line{5})/2
                 coupling_line{4}+coupling_line{3}];
        S=euler2dcm(pi*[coupling_line{6} coupling_line{7} coupling_line{8}]/180);
        inter.coupling.matrix{coupling_line{1},coupling_line{2}}=inter.coupling.matrix{coupling_line{1},coupling_line{2}}+S*diag(eigvals)*S';
        
    elseif (length(current_line)>=10)&&strcmp(current_line(1:10),'quadrupole')
        
        % Bomb out if the quadrupole tensor is not in Hz
        if nnz(current_line(11:end)=='p')~=0
            error('only Hz specification (# # # # # # #) is supported for quadrupole tensors.');
        end
        
        % Get the quadrupole tensor information
        quad_line=textscan(current_line(11:end),'%n %n %n %n %n %n %n');
        eigvals=[-quad_line{3}*(1+quad_line{4})/2
                 -quad_line{3}*(1-quad_line{4})/2
                  quad_line{3}];
        S=euler2dcm(pi*[quad_line{5} quad_line{6} quad_line{7}]/180);
        inter.coupling.matrix{quad_line{1},quad_line{1}}=inter.coupling.matrix{quad_line{1},quad_line{1}}+S*diag(eigvals)*S';
        
    end
    
end

end

% Consistency enforcement
function grumble(filename)
if ~ischar(filename)
    error('filename must be a character string.');
end
end

% According to a trade rumour, the 1983 vacancy for a Lecturer post
% at Corpus Christi College, Oxford, was contested, amongst others,
% by three candidates: Geoffrey Bodenhausen, Peter Hore and Malcolm
% Levitt. Scrolling the Scopus database to that year and sorting the
% results by citation numbers would provide much food for thought to
% any young researcher seeking to further his academic career.

