% Ellipsoid plots of chemical shielding tensors and their eigensystems.
%
% <http://spindynamics.org/wiki/index.php?title=Cst_display.m>

function cst_display(props,atoms,scaling_factor,conmatrix)

% Check consistency
grumble(props,atoms,scaling_factor,conmatrix);

% Set up graphics
figure(); hold on; colormap hot; opengl software;
light('Position',[-2,2,20]); light('Position',[10,10,10]);

% Draw the molecule
molplot(props.std_geom,conmatrix);

% Get a unit sphere
k=4; n=2^k-1; theta=pi*(-n:2:n)/n; phi=(pi/2)*(-n:2:n)'/n;
X=cos(phi)*cos(theta); Y=cos(phi)*sin(theta);
Z=sin(phi)*ones(size(theta));

% Loop over atoms
for n=1:props.natoms

    % Ignore atoms that are not on the user-supplied list
    if ismember(props.symbols{n},atoms)
        
        % Diagonalize the shielding tensor
        [eigvecs,eigvals]=eig(props.cst{n}); eigvals=diag(eigvals);
        
        % Do not plot tensors smaller than 0.1 ppm
        if norm(eigvals)>0.1
            
            % Stretch and scale sphere vertex coordinates
            coords=[reshape(X*eigvals(1)*scaling_factor,1,length(phi)*length(theta));...
                    reshape(Y*eigvals(2)*scaling_factor,1,length(phi)*length(theta));...
                    reshape(Z*eigvals(3)*scaling_factor,1,length(phi)*length(theta))];
            
            % Rotate sphere vertex coordinates
            coords=eigvecs*coords;
            
            % Fold back sphere vertex coordinates
            Xd=reshape(coords(1,:),length(phi),length(theta));
            Yd=reshape(coords(2,:),length(phi),length(theta));
            Zd=reshape(coords(3,:),length(phi),length(theta));
            
            % Color vertices by amplitude
            Cd=sqrt(Xd.^2+Yd.^2+Zd.^2);
            
            % Translate vertices
            Xd=Xd+props.std_geom(n,1);
            Yd=Yd+props.std_geom(n,2);
            Zd=Zd+props.std_geom(n,3);
            
            % Draw the surface
            surf(Xd,Yd,Zd,Cd,'FaceAlpha',0.5,'EdgeAlpha',0.1);
            
            % Scale eigenvectors
            vector_a=eigvecs(:,1)*eigvals(1)*scaling_factor;
            vector_b=eigvecs(:,2)*eigvals(2)*scaling_factor;
            vector_c=eigvecs(:,3)*eigvals(3)*scaling_factor;
            
            % Colour eigenvectors
            col_a='r-'; col_b='r-'; col_c='r-';
            if eigvals(1)<0, col_a='b-'; end
            if eigvals(2)<0, col_b='b-'; end
            if eigvals(3)<0, col_c='b-'; end
            
            % Draw eigenvectors
            plot3([props.std_geom(n,1)-vector_a(1) props.std_geom(n,1)+vector_a(1)],...
                  [props.std_geom(n,2)-vector_a(2) props.std_geom(n,2)+vector_a(2)],...
                  [props.std_geom(n,3)-vector_a(3) props.std_geom(n,3)+vector_a(3)],col_a,'LineWidth',1);
            plot3([props.std_geom(n,1)-vector_b(1) props.std_geom(n,1)+vector_b(1)],...
                  [props.std_geom(n,2)-vector_b(2) props.std_geom(n,2)+vector_b(2)],...
                  [props.std_geom(n,3)-vector_b(3) props.std_geom(n,3)+vector_b(3)],col_b,'LineWidth',1);
            plot3([props.std_geom(n,1)-vector_c(1) props.std_geom(n,1)+vector_c(1)],...
                  [props.std_geom(n,2)-vector_c(2) props.std_geom(n,2)+vector_c(2)],...
                  [props.std_geom(n,3)-vector_c(3) props.std_geom(n,3)+vector_c(3)],col_c,'LineWidth',1);
            
        end
    end
end

% Tidy up the picture
axis square; axis tight; axis equal;
set(gca,'Projection','perspective');
box on; grid on;

end

% Consistency enforcement
function grumble(props,atoms,scaling_factor,conmatrix)
if ~isfield(props,'std_geom')
    error('props structure does not contain std_geom field with atomic coordinates.');
end
if ~isfield(props,'cst')
    error('props structure does not contain cst field with chemical shielding tensors.');
end
if ~iscell(atoms)
    error('atoms parameter must be a cell array of character strings.');
end
if (~isnumeric(scaling_factor))||(~isreal(scaling_factor))
    error('scaling_factor must be a real number.');
end
if ~isempty(conmatrix)
    if ~all(size(conmatrix)==[size(std_geom,1) size(std_geom,1)])
        error('both dimensions of conmatrix must be equal to the number of atoms in the system.');
    end
end
end

% I'll be more enthusiastic about encouraging thinking outside the box when
% there's evidence of any thinking going on inside it.
%
% Terry Pratchett

