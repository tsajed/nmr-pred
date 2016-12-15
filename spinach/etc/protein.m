% Protein data import function.
%
% <http://spindynamics.org/wiki/index.php?title=Protein.m>

function [sys,inter]=protein(pdb_file,bmrb_file,options)

% Check consistency
grumble(pdb_file,bmrb_file,options);

% Parse the PDB file
[pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords]=read_pdb_pro(pdb_file,options.pdb_mol);

% Parse the BMRB file
[bmrb_aa_num,bmrb_aa_typ,bmrb_atom_id,bmrb_chemsh]=read_bmrb(bmrb_file);

% Remove oxygens, sulphurs and terminal atoms
kill_mask=ismember(pdb_atom_id,{'O','OE','OE1','OE2','OD1','OD2','OG','OG1','HG1',...
                                'OG2','OH','HH','SD','SG','OXT','O''','O'''''});
pdb_aa_num(kill_mask)=[]; pdb_atom_id(kill_mask)=[]; 
pdb_aa_typ(kill_mask)=[]; pdb_coords(kill_mask)=[];
disp('WARNING: oxygen, sulphur, and OH protons will not appear in the simulation.');

% Match chemical shifts
pdb_chemsh=cell(numel(pdb_atom_id),1);
disp(' '); disp('############# CHEMICAL SHIFT IMPORT SUMMARY ###############');
for n=1:numel(pdb_atom_id)
    
    % Pull the current amino acid from BMRB
    select_mask=(bmrb_aa_num==pdb_aa_num(n));
    bmrb_atoms=bmrb_atom_id(select_mask);
    bmrb_shifts=bmrb_chemsh(select_mask);
    
    % Make sure amino acid types match
    if ~all(strcmp(pdb_aa_typ{n},bmrb_aa_typ(select_mask)))
        disp(pdb_aa_typ{n}); disp(bmrb_aa_typ(select_mask));
        error('Amino acid type mismatch between PDB and BMRB files.');
    end
    
    % Ugly heuristics (sorry!) to match PDB and BMRB data
    if ismember(pdb_atom_id{n},bmrb_atoms)
        pdb_chemsh{n}=bmrb_shifts(strcmp(pdb_atom_id{n},bmrb_atoms));
        disp(['Imported chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif strcmp(pdb_atom_id{n},'N')&&(pdb_aa_num(n)==1)
        disp(['WARNING: N-terminal NH2 group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    elseif ismember(pdb_atom_id{n},{'HG'})&&ismember(pdb_aa_typ(n),{'SER'})
        disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    elseif ismember(pdb_atom_id{n},{'HG1'})&&ismember(pdb_aa_typ(n),{'THR'})
        disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    elseif ismember(pdb_atom_id{n},{'HB3'})&&ismember('HB2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'MET','GLU','GLN','LYS','LEU','SER','HIS','ARG'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HB2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HB1','HB2','HB3'})&&ismember('HB',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ALA'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HB',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HG21','HG22','HG23','HG1','HG3'})&&ismember('HG2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ILE','THR','VAL','GLU','LYS','PRO','GLN','ARG'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HG2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HG11','HG12','HG13','HG2','HG3'})&&ismember('HG1',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'VAL'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HG1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HG13'})&&ismember('HG12',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ILE'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HG12',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HD11','HD12','HD13','HD2','HD3'})&&ismember('HD1',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'ILE','LEU'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HD1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HD21','HD22','HD23','HD1','HD3'})&&ismember('HD2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'LYS','LEU','PRO','ARG'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HD2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HE3'})&&ismember('HE2',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'LYS'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HE2',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HE1','HE2','HE3'})&&ismember('HE',bmrb_atoms)&&ismember(pdb_aa_typ{n},{'MET'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HE',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'CD2'})&&ismember('CD1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('CD1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'CE2'})&&ismember('CE1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('CE1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HD2'})&&ismember('HD1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HD1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HE2'})&&ismember('HE1',bmrb_atoms)&&ismember(pdb_aa_typ(n),{'PHE','TYR'})
        pdb_chemsh{n}=bmrb_shifts(strcmp('HE1',bmrb_atoms));
        disp(['WARNING: replicated chemical shift for ' pad([pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ':'],15,'right') pad(num2str(pdb_chemsh{n},'%6.3f'),6,'left') ' ppm']);
    elseif ismember(pdb_atom_id{n},{'HH'})&&ismember(pdb_aa_typ(n),{'TYR'})
        disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    elseif ismember(pdb_atom_id{n},{'HD1'})&&ismember(pdb_aa_typ(n),{'HIS'})
        disp(['WARNING: ring NH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    elseif ismember(pdb_atom_id{n},{'NZ','HZ1','HZ2','HZ3'})&&ismember(pdb_aa_typ(n),{'LYS'})
        disp(['WARNING: side chain NH3+ group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    elseif ismember(pdb_atom_id{n},{'NE','CZ','NH1','NH2','HE','HH11','HH12','HH21','HH22'})&&ismember(pdb_aa_typ(n),{'ARG'})
        disp(['WARNING: side chain nitrogen block of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.']);
    end
    
end

% Estimate J-couplings
scalar_couplings=guess_j_pro(pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords);

% Estimate chemical shielding anisotropies
CSAs=guess_csa_pro(pdb_aa_num,pdb_atom_id,pdb_coords);

% Assign isotopes and labels
isotopes=cell(1,numel(pdb_atom_id));
for n=1:numel(pdb_atom_id) 
    switch pdb_atom_id{n}(1)
        case 'H'
            isotopes{n}='1H';
        case 'C'
            isotopes{n}='13C';
        case 'N'
            isotopes{n}='15N';
        otherwise
            error('unknown atom type.');
    end
end

% Process atom selection specification
if isnumeric(options.select)
    
    % Numeric selection needs no special processing

elseif strcmp(options.select,'backbone')
    
    % Import backbone up to CB and HB
    subset=ismember(pdb_atom_id,{'H','N','C','CA','HA','CB','HB','HB1','HB2','HB3'});
    
elseif strcmp(options.select,'backbone-minimal')
    
    % Import minimal backbone
    subset=ismember(pdb_atom_id,{'H','N','C','CA','HA'});
    
elseif strcmp(options.select,'backbone-hsqc')
    
    % Import backbone up to CB and HB and also ASN and GLN amide groups
    subset=ismember(pdb_atom_id,{'H','N','C','CA','HA','NE2','HE21','HE22','CD','CG',...
                                 'ND2','HD21','HD22','CB','HB','HB1','HB2','HB3'});
                             
elseif strcmp(options.select,'all')
    
    % Import everything
    subset=true(size(pdb_atom_id));

else 
    
    % Complain and bomb out
    error('incorrect subset selection specification.');
    
end

% Find missing chemical shifts
missing_shifts=find(cellfun(@isempty,pdb_chemsh))';
disp(' '); disp('############# SUMMARY OF MISSING ASSIGNMENTS ###############');

% Process unassigned chemical shifts
if strcmp(options.noshift,'keep')
    
    % Put unassigned spins between -1.0 and 0.0 ppm
    erzatz_shifts=linspace(-1,0,numel(missing_shifts));
    for n=1:numel(missing_shifts)
        disp(['Missing assignment of ' pdb_aa_typ{missing_shifts(n)} '(' num2str(pdb_aa_num(missing_shifts(n)))...
              ')-' pdb_atom_id{missing_shifts(n)} ': the spin is placed at ' num2str(erzatz_shifts(n)) ' ppm.']);
        pdb_chemsh{missing_shifts(n)}=erzatz_shifts(n);
    end
    
elseif strcmp(options.noshift,'delete')
    
    % Delete unassigned spins
    for n=missing_shifts
        disp(['Missing assignment of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')-' pdb_atom_id{n} ...
              ': the atom will not appear in the simulation.']);
    end
    subset=subset&(~cellfun(@isempty,pdb_chemsh));
    
else
    
    % Complain and bomb out
    error('incorrect value of options.noshift parameter.');
    
end

% Apply the selection
sys.isotopes=isotopes(subset);
sys.labels=pdb_atom_id(subset)';
inter.zeeman.scalar=pdb_chemsh(subset)';
inter.zeeman.matrix=CSAs(subset)';
inter.coupling.scalar=scalar_couplings(subset,subset);
inter.coordinates=pdb_coords(subset);

end

% Consistency enforcement
function grumble(pdb_file,bmrb_file,options)
if ~ischar(pdb_file)
    error('pdb_file must be a character string specifying a file name.');
end
if ~ischar(bmrb_file)
    error('bmrb_file must be a character string specifying a file name.');
end
if ~isfield(options,'select')
    error('options.select switch must be specfied.');
elseif (~isnumeric(options.select))&&(~ismember(options.select,{'all','backbone','backbone-hsqc','backbone-minimal'}))
    error('invalid value for options.select, please refer to the manual.');
end
if ~isfield(options,'pdb_mol')
    error('options.pdb_mol switch must be specfied.');
elseif (~isnumeric(options.pdb_mol))
    error('invalid value for options.pdb_mol, please refer to the manual.');
end
if ~isfield(options,'noshift')
    error('options.noshift switch must be specfied.');
elseif (~isnumeric(options.noshift))&&(~ismember(options.noshift,{'keep','delete'}))
    error('invalid value for options.noshift, please refer to the manual.');
end
end

% We saw that we'd been given a law to live by, a moral law, they called
% it, which punished those who observed it -- for observing it. The more 
% you tried to live up to it, the more you suffered; the more you cheated
% it, the bigger reward you got.
%
% Ayn Rand, "Atlas Shrugged"

