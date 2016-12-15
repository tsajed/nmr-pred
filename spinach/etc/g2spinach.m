% Forms Spinach data structures from Gaussian parsing output returned by gparse.m function.
%
% <http://spindynamics.org/wiki/index.php?title=G2spinach.m>

function [sys,inter]=g2spinach(props,particles,references,options)

% Check consistency
grumble(props,particles,references,options);

% Index the particles to include
sys.isotopes={}; index=[]; ref_index=[];
for n=1:length(props.symbols) %#ok<*AGROW>
    for k=1:length(particles) 
        if strcmp(props.symbols{n},particles{k}{1})
            sys.isotopes=[sys.isotopes particles{k}{2}];
            index=[index n]; ref_index=[ref_index k];
        end
    end
end
nspins=length(index);

% Decide whether to include coordinates
if nargin<4
    include_xyz=true();
else
    if ~isfield(options,'no_xyz')
        include_xyz=true();
    elseif ~options.no_xyz
        include_xyz=true();
    else
        include_xyz=false();
    end
end

% Process coordinates
if include_xyz
    inter.coordinates=props.std_geom(index,:);
    inter.coordinates=num2cell(inter.coordinates,2);
end

% Decide which magnetic parameters to return
switch ismember('E',[particles{:}])
    
    % EPR parameterization
    case 1
        
        % Add the electron as the last spin in the isotope list
        sys.isotopes=[sys.isotopes, 'E']; nspins=nspins+1;
        
        % Electron coordinates should be treated as unknown
        if include_xyz, inter.coordinates{end+1}=[]; end
        
        % All Zeeman tensors are zero except for the g-tensor of the electron
        inter.zeeman.matrix=cell(1,nspins);
        inter.zeeman.matrix{nspins}=props.g_tensor.matrix;
        
        % All couplings are zero except for the hyperfine couplings to the electron
        inter.coupling.matrix=mat2cell(zeros(3*nspins,3*nspins),3*ones(nspins,1),3*ones(nspins,1));
        for n=1:(nspins-1)
            inter.coupling.matrix{n,end}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
            inter.coupling.matrix{end,n}=1e6*gauss2mhz(props.hfc.full.matrix{index(n)}/2);
        end
        
        % Remove small hyperfine couplings
        if (nargin==4)&&isfield(options,'min_hfc')
            for n=1:nspins
                for k=1:nspins
                    if norm(inter.coupling.matrix{n,k},'fro')<options.min_hfc
                        inter.coupling.matrix{n,k}=[];
                    end
                end
            end
        end
        
        % Remove the nuclei that are not coupled to the electron
        if (nargin==4)&&isfield(options,'purge')&&strcmp(options.purge,'on')
            killing_pattern=cellfun(@isempty,inter.coupling.matrix(:,nspins)); killing_pattern(end)=0;
            inter.coupling.matrix(:,killing_pattern)=[];
            inter.coupling.matrix(killing_pattern,:)=[];
            inter.zeeman.matrix(killing_pattern)=[];
            sys.isotopes(killing_pattern)=[];
        end
        
    % NMR parameterization    
    case 0
        
        % Reference to bare nuclei if the user did not supply any reference values
        if (nargin<3)||(~any(references)), references=zeros(size(particles)); end
        
        % Absorb Zeeman tensors
        if isfield(props,'cst')
            inter.zeeman.matrix=cell(1,nspins);
            for n=1:nspins
                inter.zeeman.matrix{n}=-props.cst{index(n)}+eye(3)*references(ref_index(n));
            end
        end
        
        % Absorb quadrupolar couplings
        if isfield(props,'nqi')
            inter.coupling.matrix=cell(nspins,nspins);
            for n=1:nspins
                inter.coupling.matrix{n,n}=props.nqi{index(n)};
            end
        end
        
        % Absorb spin-rotation couplings
        if isfield(props,'src')
            inter.spinrot.matrix=cell(nspins);
            for n=1:nspins
                inter.spinrot.matrix{n}=props.src{index(n)};
            end
        end
        
        % Absorb and prune scalar couplings
        if isfield(props,'j_couplings')
            inter.coupling.scalar=props.j_couplings(index,index)/2;
            if (nargin==4)&&isfield(options,'min_j')
                inter.coupling.scalar=inter.coupling.scalar.*(abs(inter.coupling.scalar)>options.min_j);
            end
            inter.coupling.scalar=num2cell(inter.coupling.scalar);
        end
       
end

end

% Consistency enforcement
function grumble(props,nuclei,references,options) %#ok<INUSD>
if ~isstruct(props)
    error('the first argument must be a structure returned by gparse().');
end
if ~iscell(nuclei)
    error('the second argument must have the form {{''H'',''1H''},{''C'',''13C''},...}');
end
if (~isnumeric(references))||(numel(nuclei)~=numel(references))
    error('references must be a numerical array with the same number of entries as nuclei.');
end
end

% ACHTUNG! ALLES LOOKENSPEEPERS
% Das Computermachine ist nicht fur gefingerpoken und mittengrabben.
% Ist easy schnappen der Springenwerk, blowenfusen und poppencorken
% mit Spitzensparken. Ist nicht fur gewerken bei das Dumpkopfen. Die
% rubbernecken Sichtseeren keepen Hands in die Pockets muss, relaxen
% und watchen die Blinkenlichten.

