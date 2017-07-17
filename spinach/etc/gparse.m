% A parser for Gaussian03 and Gaussian09 calculation logs.
%
% <http://spindynamics.org/wiki/index.php?title=Gparse.m>

function props=gparse(filename,options)

% Read the file
file_id=fopen(filename,'r');
g03_output=textscan(file_id,'%s','delimiter','\n');
fclose(file_id); g03_output=g03_output{1};

% Deblank all lines
for n=1:numel(g03_output), g03_output(n)=deblank(g03_output(n)); end

% Refuse to process files without #p option
for n=1:numel(g03_output)
    current_line=char(g03_output(n));
    if (numel(current_line)>=2)&&strcmp(current_line(1:2),'# ')
        error('This Gaussian log has the detailed printing option switched off. Re-run with #p in the route line.');
    end
end

% Set the default error flag
props.error=0;

% Set the default complete flag
props.complete=0;

% Parse the file
for n=1:length(g03_output)
   
   % Read the input orientation and atomic numbers
   if strcmp(g03_output(n),'Input orientation:')
      k=n+5; m=1;
      while ~strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')
         S=eval(['[' char(g03_output(k)) ']']); atoms(m,:)=S(4:6); atomic_numbers(m)=S(2); k=k+1; m=m+1; %#ok<AGROW>
      end
      natoms=m-1; props.inp_geom=atoms; props.natoms=natoms;
      disp('Gaussian import: found input orientation.');
   end

   % Read the standard orientation and atomic numbers
   if strcmp(g03_output(n),'Standard orientation:')
      k=n+5; m=1;
      while ~strcmp(deblank(g03_output(k)),'---------------------------------------------------------------------')
         S=eval(['[' char(g03_output(k)) ']']); atoms(m,:)=S(4:6); atomic_numbers(m)=S(2); k=k+1; m=m+1;
      end
      natoms=m-1; props.std_geom=atoms; props.natoms=natoms;
      disp('Gaussian import: found standard orientation.');
   end

   % Read the SCF energy
   current_line=char(g03_output(n));
   if length(current_line)>10 && strcmp(current_line(1:9),'SCF Done:')
       scan_data=textscan(current_line,'SCF Done: %s = %f','Delimiter',' ','MultipleDelimsAsOne',1);
       props.method=char(scan_data{1}); props.method=props.method(3:(end-1)); props.energy=scan_data{2};
       disp('Gaussian import: found SCF energy.');
   end

   % Read isotropic hyperfine couplings
   if strcmp(g03_output(n),'Isotropic Fermi Contact Couplings')
      props.hfc.iso=zeros(natoms,1);
      for k=1:natoms
         S=char(g03_output(n+1+k)); S=eval(['[' S(20:end) ']']); props.hfc.iso(k)=S(3);
      end
      disp('Gaussian import: found isotropic hyperfine couplings.');
   end

   % Read and symmetrize anisotropic hyperfine couplings
   if strcmp(g03_output(n),'Anisotropic Spin Dipole Couplings in Principal Axis System')
      props.hfc.full.eigvals=cell(natoms,1);
      props.hfc.full.eigvecs=cell(natoms,1);
      props.hfc.full.matrix=cell(natoms,1);
      for k=1:natoms
         baa=g03_output(n+4*k+1); baa=char(baa); baa=eval(['[' baa(4:end) ']']);
         bbb=g03_output(n+4*k+2); bbb=char(bbb); bbb=eval(['[' bbb(15:end) ']']);
         bcc=g03_output(n+4*k+3); bcc=char(bcc); bcc=eval(['[' bcc(4:end) ']']);
         props.hfc.full.eigvals{k}=[baa(3) bbb(3) bcc(3)]+props.hfc.iso(k);
         props.hfc.full.eigvecs{k}=[baa(5:7)' bbb(5:7)' bcc(5:7)'];
         props.hfc.full.matrix{k}=props.hfc.full.eigvecs{k}*diag(props.hfc.full.eigvals{k})*props.hfc.full.eigvecs{k}';
         if (~exist('options','var'))||(~ismember('hfc_nosymm',options))
            props.hfc.full.matrix{k}=(props.hfc.full.matrix{k}+props.hfc.full.matrix{k}')/2;
         end
      end
      disp('Gaussian import: found anisotropic hyperfine couplings.');
   end

   % Read and symmetrize the g-tensor (Gaussian03 for Unix)
   if strcmp(g03_output(n),'g tensor [g = g_e + g_RMC + g_DC + g_OZ/SOC]:')
      line1=char(deblank(g03_output(n+1))); line2=char(deblank(g03_output(n+2))); line3=char(deblank(g03_output(n+3)));
      g=eval(['[' line1(5:18) '   ' line1(26:39) '   ' line1(46:60) ';   '...
                  line2(5:18) '   ' line2(26:39) '   ' line2(46:60) ';   '...
                  line3(5:18) '   ' line3(26:39) '   ' line3(46:60) ']']);
      if (~exist('options','var'))||(~ismember('g_nosymm',options))
         g=(g+g')/2;
      end
      [V,D]=eig(g);
      props.g_tensor.eigvecs=V;
      props.g_tensor.eigvals=diag(D)';
      props.g_tensor.matrix=g;
      disp('Gaussian import: found g-tensor.');
   end
   
   % Read and symmetrize the g-tensor (Gaussian03 for Windows)
   if strcmp(g03_output(n),'g tensor (ppm):')
      line1=char(deblank(g03_output(n+1))); line2=char(deblank(g03_output(n+2))); line3=char(deblank(g03_output(n+3)));
      g=eval(['[' line1(4:15) '   ' line1(23:35) '   ' line1(42:54) ';   '...
                  line2(4:15) '   ' line2(23:35) '   ' line2(42:54) ';   '...
                  line3(4:15) '   ' line3(23:35) '   ' line3(42:54) ']']);
      if (~exist('options','var'))||(~ismember('g_nosymm',options))
         g=(g+g')/2;
      end
      [V,D]=eig(g);
      props.g_tensor.eigvecs=V;
      props.g_tensor.eigvals=diag(D)';
      props.g_tensor.matrix=g;
      disp('Gaussian import: found g-tensor.');
   end
   
   % Read and symmetrize chemical shielding tensors
   if strcmp(g03_output(n),'SCF GIAO Magnetic shielding tensor (ppm):')||strcmp(g03_output(n),'Magnetic shielding (ppm):')
       cst=cell(natoms,1); 
       for k=1:natoms
           line1=char(g03_output(n+5*k-3)); line2=char(g03_output(n+5*k-2)); line3=char(g03_output(n+5*k-1));
           cst{k}=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                   eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                   eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])];
           if (~exist('options','var'))||(~ismember('cst_nosymm',options))
               cst{k}=(cst{k}+cst{k}')/2;
           end
       end
       props.cst=cst;
       disp('Gaussian import: found chemical shielding tensors.');
   end
   
   % Read scalar couplings
   if strcmp(g03_output(n),'Total nuclear spin-spin coupling J (Hz):')
       j_couplings=zeros(natoms,ceil(natoms/5)*5);
       for k=1:ceil(natoms/5)
           lines_to_read=natoms-(k-1)*5;
           current_block=zeros(natoms,5);
           for m=(n+2):(n+1+lines_to_read)
               current_line=eval(['[' char(g03_output(m)) ']']);
               current_block(current_line(1),1:(length(current_line)-1))=current_line(2:end);
           end
           j_couplings(:,(5*(k-1)+1):(5*k))=current_block;
           n=n+lines_to_read+1; %#ok<FXSET>
       end
       j_couplings=j_couplings(1:natoms,1:natoms); j_couplings=j_couplings+j_couplings';
       props.j_couplings=j_couplings;
       disp('Gaussian import: found isotropic J-couplings.');
   end
   
   % Read spin-rotation couplings
   if strcmp(g03_output(n),'nuclear spin - molecular rotation tensor [C] (MHz):')
       props.srt=cell(natoms,1);
       for k=1:natoms
           line3=char(g03_output(n+4*k));   line2=char(g03_output(n+4*k-1));
           line1=char(g03_output(n+4*k-2)); line0=char(g03_output(n+4*k-3));
           atom_num=regexp(line0,'^[0-9]*','match'); atom_num=eval(atom_num{1});
           props.srt{atom_num}=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                                eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                                eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])]*1e6;
           if strcmp(g03_output(n+4*k+1),'Dipole moment (Debye):'), break; end
           if strcmp(g03_output(n+4*k+1),'Nuclear quadrupole coupling constants [Chi] (MHz):'), break; end
       end
       disp('Gaussian import: found nuclear spin-rotation tensors.');
   end
   
   % Read quadrupole couplings and kill their trace
   if strcmp(g03_output(n),'Nuclear quadrupole coupling constants [Chi] (MHz):')
       props.nqi=cell(natoms,1);
       for k=1:natoms
           line3=char(g03_output(n+4*k));   line2=char(g03_output(n+4*k-1)); 
           line1=char(g03_output(n+4*k-2)); line0=char(g03_output(n+4*k-3));
           atom_num=regexp(line0,'^[0-9]*','match'); atom_num=eval(atom_num{1});
           props.nqi{atom_num}=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                                eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                                eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])]*1e6;
           props.nqi{atom_num}=props.nqi{atom_num}-eye(3)*trace(props.nqi{atom_num})/3;
           if strcmp(g03_output(n+4*k+1),'Dipole moment (Debye):'), break; end
       end
       disp('Gaussian import: found nuclear quadrupole tensors.');
   end
   
   % Read magnetic susceptibility tensor (GIAO)
   if strcmp(g03_output(n),'Magnetic susceptibility tensor (cgs-ppm):')
       line1=char(g03_output(n+1)); line2=char(g03_output(n+2)); line3=char(g03_output(n+3));
       props.chi=[eval(['[' line1(4:19) '     ' line1(25:40) '     ' line1(46:end) ']']);
                  eval(['[' line2(4:19) '     ' line2(25:40) '     ' line2(46:end) ']']);
                  eval(['[' line3(4:19) '     ' line3(25:40) '     ' line3(46:end) ']'])];
       props.chi=cgsppm2ang(props.chi);
       disp('Gaussian import: found magnetic susceptibility.');
   end
   
   % Read magnetic susceptibility tensor (CSGT)
   if strcmp(g03_output(n),'Magnetic susceptibility (cgs-ppm):')
       line1=char(g03_output(n+2)); line2=char(g03_output(n+3)); line3=char(g03_output(n+4));
       props.chi=[eval(['[' line1(4:14) '     ' line1(21:31) '     ' line1(38:48) ']']);
                  eval(['[' line2(4:14) '     ' line2(21:31) '     ' line2(38:48) ']']);
                  eval(['[' line3(4:14) '     ' line3(21:31) '     ' line3(38:48) ']'])];
       props.chi=cgsppm2ang(props.chi);
       disp('Gaussian import: found magnetic susceptibility.');
   end
   
   % Read Gibbs free energy
   if (numel(g03_output{n})>43)&&strcmp(g03_output{n}(1:43),'Sum of electronic and thermal Free Energies')
       props.gibbs=eval(g03_output{n}(45:end));
       disp('Gaussian import: found Gibbs free energy.');
   end
   
   % Check for error flags
   if (numel(g03_output{n})>17)&&strcmp(g03_output{n}(1:17),'Error termination')
       props.error=1; disp('Gaussian import: error message detected.');
   end
   
   % Check for incomplete calculations
   if (numel(g03_output{n})>18)&&strcmp(g03_output{n}(1:18),'Normal termination')
       props.complete=1;
   end
   if (numel(g03_output{n})>24)&&strcmp(g03_output{n}(1:24),'Entering Gaussian System')
       props.complete=0;
   end
   
end
   
% Assign atomic symbols
periodic_table={'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar',...
                'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
                'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
                'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
                'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
                'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',...
                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg'};

if exist('atomic_numbers','var')
    props.symbols=periodic_table(atomic_numbers);
    props.atomic_numbers=atomic_numbers;
end

% Assign source information
props.filename=filename;

end

% Moral indignation is jealousy with a halo.
% 
% H.G. Wells
